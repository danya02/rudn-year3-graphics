use bresenham::Bresenham;
use nalgebra::{Matrix, Matrix4, OPoint, Point3, UnitQuaternion, Vector3, Vector4};

use crate::math::{
    edge_function, is_pixel_overlapping, quat_to_rotmat, sort_points_clockwise, Quat, Vec3,
};

/// Represents the 3D world to draw.
pub struct World {
    pub camera: Camera,
    pub shapes: Vec<Shape>,
}

pub struct Camera {
    pub location: Vec3,
    pub rotation: Quat,
    pub projection: Projection,
    pub fov_radians: f64,
    pub ortho_bounds: (f64, f64, f64, f64),
}

impl Camera {
    #[rustfmt::skip]
    fn world_to_camera_basic(&self, near: f64, far: f64) -> Matrix4<f64> {
        match self.projection {
            Projection::Orthographic => {
                let scale = 1.0 / f64::tan(self.fov_radians);
                let mut m: Matrix4<f64> = Matrix4::identity();
                m[(0, 0)] = scale;
                m[(1, 1)] = scale;
                m[(2, 2)] = -far / (far - near);
                m[(3, 2)] = -far * near / (far - near);
                m[(2, 3)] = -1.0;
                m[(3, 3)] = 0.0;
                m
            }
            Projection::Perspective => {
                let (t, b, l, r) = self.ortho_bounds;
                let (f, n) = (far, near);

                let m = Matrix4::new(
                        2.0/(r-l), 0.0, 0.0, 0.0,
                        0.0, 2.0/(t-b), 0.0, 0.0,
                        0.0, 0.0, -2.0/(f-n), 0.0,
                        -(r+l)/(r-l), (-t+b)/(t-b), -(f+n)/(f-n), 1.0
                    );
                m
            }
        }
    }

    pub fn world_to_camera(&self, near: f64, far: f64) -> Matrix4<f64> {
        let mut m = self.world_to_camera_basic(near, far); //* quat_to_rotmat(self.rotation);
        m[(3, 0)] = self.location[0];
        m[(3, 1)] = self.location[1];
        m[(3, 2)] = self.location[2];
        m
    }
}

pub enum Projection {
    Perspective,
    Orthographic,
}

pub struct Shape {
    pub origin: Vec3,
    pub rotation: Quat,
    pub vertices: Vec<Vec3>,
    pub triangles: Vec<(usize, usize, usize)>,
}

impl World {
    pub fn rasterize(&self, frame: &mut [u8], width: u32, height: u32) {
        for (i, pixel) in frame.chunks_exact_mut(4).enumerate() {
            let x = (i % width as usize) as usize;
            let y = (i / width as usize) as usize;

            let rgba = { [0, 0, 0, 0xff] };

            pixel.copy_from_slice(&rgba);
        }

        let mut depth_buffer = vec![vec![f64::INFINITY; width as usize]; height as usize];

        let mut set_at = |x, y, color: [u8; 4]| {
            frame
                .chunks_exact_mut(4)
                .skip(y * width as usize)
                .skip(x)
                .next()
                .unwrap()
                .copy_from_slice(&color)
        };

        let camera_matrix = self.camera.world_to_camera(0.01, 1000.0);

        println!("{camera_matrix:?}");

        for shape in self.shapes.iter() {
            for (a, b, c) in shape.triangles.iter() {
                // Transform the triangle according to the shape: first shift them, then rotate the points.
                let a = shape.vertices[*a];
                let b = shape.vertices[*b];
                let c = shape.vertices[*c];
                let a = shape.rotation.transform_point(&Point3::new(a.x, a.y, a.z));
                let b = shape.rotation.transform_point(&Point3::new(b.x, b.y, b.z));
                let c = shape.rotation.transform_point(&Point3::new(c.x, c.y, c.z));

                let a = a + shape.origin;
                let b = b + shape.origin;
                let c = c + shape.origin;

                fn promote(src: Point3<f64>) -> Vector4<f64> {
                    Vector4::new(src[0], src[1], src[2], 0.0)
                }
                fn demote(src: Vector4<f64>) -> Vec3 {
                    Vec3::new(src[0], src[1], src[2])
                }

                let mut a = promote(a);
                let mut b = promote(b);
                let mut c = promote(c);

                a.x += camera_matrix[(3, 0)];
                a.y += camera_matrix[(3, 1)];
                a.z += camera_matrix[(3, 2)];
                b.x += camera_matrix[(3, 0)];
                b.y += camera_matrix[(3, 1)];
                b.z += camera_matrix[(3, 2)];
                c.x += camera_matrix[(3, 0)];
                c.y += camera_matrix[(3, 1)];
                c.z += camera_matrix[(3, 2)];

                let mut a = camera_matrix * a;
                let mut b = camera_matrix * b;
                let mut c = camera_matrix * c;

                // println!(
                //     "{} {} {}",
                //     camera_matrix[(3, 0)],
                //     camera_matrix[(3, 1)],
                //     camera_matrix[(3, 2)]
                // );

                let a = demote(a);
                let b = demote(b);
                let c = demote(c);

                // Ensure winding direction
                let (a, b, c) = sort_points_clockwise((a, b, c));

                // Compute bounding box, clipped inside canvas
                let x0 = a[0].min(b[0].min(c[0]));
                let y0 = a[1].min(b[1].min(c[1]));
                let x1 = a[0].max(b[0].max(c[0]));
                let y1 = a[1].max(b[1].max(c[1]));

                let x0 = x0 as isize;
                let y0 = y0 as isize;
                let x1 = x1 as isize;
                let y1 = y1 as isize;

                // If the bounding box is entirely contained outside the render area, skip it.
                if x1 < 0 || x0 > width as isize || y1 < 0 || y0 > height as isize {
                    continue;
                }

                let x0 = x0.clamp(0, width as isize - 1);
                let x1 = x1.clamp(0, width as isize - 1);
                let y0 = y0.clamp(0, height as isize - 1);
                let y1 = y1.clamp(0, height as isize - 1);

                // Draw the bounding box
                // let ba = (x0, y0);
                // let bb = (x0, y1);
                // let bc = (x1, y1);
                // let bd = (x1, y0);
                // for (x, y) in Bresenham::new(ba, bb)
                //     .chain(Bresenham::new(bb, bc))
                //     .chain(Bresenham::new(bc, bd))
                //     .chain(Bresenham::new(bd, ba))
                // {
                //     set_at(x as usize, y as usize, [255, 0, 0, 255]);
                // }

                // Draw the wireframe.
                for (x, y) in Bresenham::new(
                    (a[0] as isize, a[1] as isize),
                    (b[0] as isize, b[1] as isize),
                )
                .chain(Bresenham::new(
                    (b[0] as isize, b[1] as isize),
                    (c[0] as isize, c[1] as isize),
                ))
                .chain(Bresenham::new(
                    (c[0] as isize, c[1] as isize),
                    (a[0] as isize, a[1] as isize),
                )) {
                    if x < 0 || x >= width as isize || y < 0 || y >= height as isize {
                        continue;
                    }

                    set_at(x as usize, y as usize, [255, 255, 255, 255]);
                }
                // For each pixel in the bounding box, check if it is contained inside the triangle.
                for x in x0 as usize..x1 as usize {
                    for y in y0 as usize..y1 as usize {
                        // Compute barycentric coordinates
                        let p = Vec3::new(x as f64, y as f64, 0.0);
                        let w0 = edge_function(b, c, p);
                        let w1 = edge_function(c, a, p);
                        let w2 = edge_function(a, b, p);
                        if is_pixel_overlapping((a, b, c), (w0, w1, w2)) {
                            // Interpolate Z coordinate
                            let z = 1.0 / (1.0 / a[2] * w0 + 1.0 / b[2] * w1 + 1.0 / c[2] * w2);
                            if z < depth_buffer[y][x] {
                                depth_buffer[y][x] = z;
                                set_at(x, y, [0, 0, (z * 1000.0) as u8, 255])
                            }
                        }
                    }
                }
            }
        }
    }

    pub fn update(&mut self) {
        // self.camera.location += Vec3::new(0.01, -0.01, 0.01);
        self.shapes[0].rotation =
            self.shapes[0].rotation * Quat::from_axis_angle(&Vec3::y_axis(), 0.01);
    }
}
