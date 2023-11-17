use bresenham::Bresenham;
use nalgebra::{
    Isometry3, Matrix, Matrix4, OPoint, Perspective3, Point3, Rotation3, Translation3,
    UnitQuaternion, Vector2, Vector3, Vector4,
};

use crate::math::{
    edge_function, is_pixel_overlapping, quat_to_rotmat, sort_points_clockwise,
    transform_point_by_projection, Pt3, Quat, Vec3,
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
    pub screen_size: Vector2<f64>,
}

impl Camera {
    #[rustfmt::skip]
    fn world_to_camera_basic(&self, near: f64, far: f64) -> Matrix4<f64> {
        match self.projection {
            Projection::Perspective => {
                let m = self.rotation.to_rotation_matrix().into_inner();
                let view_matrix = Matrix4::new(
                    m.m11, m.m12, m.m13, 0.0,
                    m.m21, m.m22, m.m23, 0.0,
                    m.m31, m.m32, m.m33, 0.0,
                    self.location.x, self.location.y, self.location.z, 1.0,

                );
                let projection_matrix = Perspective3::new(self.screen_size.x / self.screen_size.x, self.fov_radians, near, far);

                projection_matrix.as_matrix() * view_matrix
            }
            Projection::Orthographic => {
                let (t, b, l, r): (f64,f64,f64,f64) = todo!();
                let (f, n) = (far, near);

                let m = Matrix4::new(
                        2.0/(r-l), 0.0, 0.0, 0.0,
                        0.0, 2.0/(t-b), 0.0, 0.0,
                        (r+l)/(r-l), (t+b)/(t-b), -(f+n)/(f-n), -1.0,
                        0.0, 0.0, -(2.0*f*n)/(f-n), 1.0
                    );
                m
            }
        }
    }

    pub fn world_to_camera(&self, near: f64, far: f64) -> Matrix4<f64> {
        let mut m = self.world_to_camera_basic(near, far) * quat_to_rotmat(self.rotation);
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

        let mut set_at = |x: usize, y: usize, color: [u8; 4]| {
            if !(0..width as usize).contains(&x) {
                return;
            }
            if !(0..height as usize).contains(&y) {
                return;
            }
            frame
                .chunks_exact_mut(4)
                .skip(y * width as usize)
                .skip(x)
                .next()
                .unwrap()
                .copy_from_slice(&color)
        };

        let camera_matrix = self.camera.world_to_camera(0.01, 100.0);

        let camera_location = self.camera.location;
        let camera_direction = Vec3::x();
        let camera_direction = self.camera.rotation.transform_vector(&camera_direction);

        // println!("{camera_matrix:?}");

        let near = 0.01;
        let far = 100.0;

        for shape in self.shapes.iter() {
            let model_matrix =
                Isometry3::from_parts(Translation3::from(shape.origin), shape.rotation);
            for (a, b, c) in shape.triangles.iter() {
                // Transform the triangle according to the shape: first shift them, then rotate the points.
                let a = shape.vertices[*a];
                let b = shape.vertices[*b];
                let c = shape.vertices[*c];
                let a = Pt3::from(a);
                let b = Pt3::from(b);
                let c = Pt3::from(c);
                let a = model_matrix.transform_point(&a);
                let b = model_matrix.transform_point(&b);
                let c = model_matrix.transform_point(&c);

                // Points are now in world space

                // let a = camera_matrix.transform_point(&a);
                // let b = camera_matrix.transform_point(&b);
                // let c = camera_matrix.transform_point(&c);

                let a = transform_point_by_projection(a.coords, camera_matrix);
                let b = transform_point_by_projection(b.coords, camera_matrix);
                let c = transform_point_by_projection(c.coords, camera_matrix);

                if a.z < near || a.z > far || b.z < near || b.z > far || c.z < near || c.z > far {
                    continue;
                }

                // Points are now in camera space

                set_at(a.x as usize, a.y as usize, [255, 255, 255, 255]);
                set_at(b.x as usize, b.y as usize, [255, 255, 255, 255]);
                set_at(c.x as usize, c.y as usize, [255, 255, 255, 255]);

                // Ensure winding direction
                let (a, b, c) = sort_points_clockwise((a.into(), b.into(), c.into()));

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
                // // For each pixel in the bounding box, check if it is contained inside the triangle.
                // for x in x0 as usize..x1 as usize {
                //     for y in y0 as usize..y1 as usize {
                //         // Compute barycentric coordinates
                //         let p = Vec3::new(x as f64, y as f64, 0.0);
                //         let w0 = edge_function(b, c, p);
                //         let w1 = edge_function(c, a, p);
                //         let w2 = edge_function(a, b, p);
                //         if is_pixel_overlapping((a, b, c), (w0, w1, w2)) {
                //             // Interpolate Z coordinate
                //             let z = 1.0 / (1.0 / a[2] * w0 + 1.0 / b[2] * w1 + 1.0 / c[2] * w2);
                //             if z < depth_buffer[y][x] {
                //                 depth_buffer[y][x] = z;
                //                 set_at(x, y, [0, 0, (z * 1000.0) as u8, 255])
                //             }
                //         }
                //     }
                // }
            }
        }
    }

    pub fn update(&mut self) {
        // self.camera.location += Vec3::new(0.01, -0.01, 0.01);
        self.shapes[0].rotation =
            self.shapes[0].rotation * Quat::from_axis_angle(&Vec3::y_axis(), 0.01);
    }
}
