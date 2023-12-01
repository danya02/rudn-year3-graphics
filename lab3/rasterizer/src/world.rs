use std::convert::Infallible;

use bresenham::Bresenham;
use embedded_graphics::{
    draw_target::DrawTarget,
    geometry::{OriginDimensions, Point},
    mono_font::{
        ascii::{FONT_10X20, FONT_6X10},
        MonoTextStyle,
    },
    pixelcolor::Rgb888,
    primitives::{PointsIter, Primitive, PrimitiveStyle, PrimitiveStyleBuilder, Triangle},
    Drawable, Pixel,
};
use nalgebra::{
    Isometry3, Matrix, Matrix4, OPoint, Orthographic3, Perspective3, Point3, Rotation3, Scale3,
    Translation3, UnitQuaternion, Vector2, Vector3, Vector4,
};
use rand::{Rng, SeedableRng};

use crate::math::{edge_function, is_pixel_overlapping, sort_points_ccw, Mat4, Pt3, Quat, Vec3};

/// Represents the 3D world to draw.
pub struct World {
    pub camera: Camera,
    pub shapes: Vec<Shape>,
    pub world_bbox: Option<(Pt3, Pt3)>,
}

pub struct Camera {
    pub location: Vec3,
    pub rotation: Quat,
    pub rotation_angles: Vec3,
    pub projection: Projection,
    pub fov_radians: f64,
    pub target_height: f64,
}

impl Camera {
    #[rustfmt::skip]
    fn ortho_raw(&self, l: f64, r: f64, t: f64, b: f64, n: f64, f: f64) -> Matrix4<f64> {

        let o = 0.0;
        // Matrix4::new(
        //     2.0/(r-l), o, o, o,
        //     o, 2.0/(t-b), o, o,
        //     o, o, -2.0/(f-n), o,
        //     -(r+l)/(r-l), -(t+b)/(t-b), -(f+n)/(f-n), 1.0
        // )

        Orthographic3::new(l,r,b,t,n,f).into()
    }

    fn world_to_camera(&self) -> Mat4 {
        // Camera is at point X. Camera moves to origin -> point moves to -X
        (Rotation3::from(self.rotation.inverse()) * Translation3::from(-self.location)).into()
    }

    fn perspective_raw(&self, aspect: f64, near: f64, far: f64) -> Mat4 {
        Perspective3::new(aspect, self.fov_radians, near, far).into()
    }
}

pub enum Projection {
    Perspective,
    Orthographic,
}

pub struct Shape {
    pub origin: Pt3,
    pub rotation: Quat,
    pub scale_factor: f64,
    pub vertices: Vec<Pt3>,
    pub triangles: Vec<[usize; 3]>,
    pub style: RenderingStyle,
}

pub enum RenderingStyle {
    /// Draw the lines corresponding to the shape's edges.
    ///
    /// This actually draws the entire triangle,
    /// and so will fill the triangle if that is supplied.
    /// However, filling the triangle will result in faces being shown in an unnatural ordering.
    Wireframe { line_style: PrimitiveStyle<Rgb888> },

    /// Draw each face in a random color, and use the Z-buffer for rendering.
    UnlitRandom,
}

impl Shape {
    pub fn shape_to_world(&self) -> Mat4 {
        Mat4::from(Translation3::from(self.origin))
            * Mat4::from(Rotation3::from(self.rotation))
            * Mat4::from(Scale3::new(
                self.scale_factor,
                self.scale_factor,
                self.scale_factor,
            ))
    }

    pub fn with_location(mut self, loc: Pt3) -> Self {
        self.origin = loc;
        self
    }
}

impl World {
    pub fn rasterize<D: DrawTarget<Color = Rgb888, Error = Infallible> + OriginDimensions>(
        &mut self,
        frame: &mut D,
    ) -> Result<(), Infallible> {
        frame.clear(Rgb888::new(0, 0, 0))?;

        let mut rng = rand::rngs::SmallRng::seed_from_u64(0);

        let (mut maxx, mut maxy, mut maxz) =
            (f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        let (mut minx, mut miny, mut minz) = (f64::INFINITY, f64::INFINITY, f64::INFINITY);

        let min;
        let max;

        if let Some((minb, maxb)) = self.world_bbox {
            min = minb;
            max = maxb;
        } else {
            for shape in self.shapes.iter() {
                let shape_to_world = shape.shape_to_world();
                for tri in shape.triangles.iter() {
                    for point in tri.iter().map(|idx| shape.vertices[*idx]) {
                        let point = shape_to_world.transform_point(&point);
                        maxx = maxx.max(point.x);
                        maxy = maxy.max(point.y);
                        maxz = maxz.max(point.z);

                        minx = minx.min(point.x);
                        miny = miny.min(point.y);
                        minz = minz.min(point.z);
                    }
                }
            }
            min = Pt3::new(minx, miny, minz);
            max = Pt3::new(maxx, maxy, maxz);
            self.world_bbox = Some((min, max));
        }
        // let min = w.transform_point(&min);
        // let max = w.transform_point(&max);

        let (minx, miny, minz) = (min.x, min.y, min.z);
        let (maxx, maxy, maxz) = (max.x, max.y, max.z);

        let max = [minx.abs(), maxx.abs(), miny.abs(), maxy.abs()]
            .into_iter()
            .fold(f64::NEG_INFINITY, f64::max);
        let width = frame.size().width;
        let height = frame.size().height;
        let fwidth = width as f64;
        let fheight = height as f64;
        let iwidth = width as i32;
        let iheight = height as i32;

        let aspect_ratio = fwidth / fheight;
        let r = max * aspect_ratio;
        let l = -r;
        let t = max;
        let b = -t;

        // let camera_to_screen = self.camera.ortho_raw(r, l, t, b, 0.01, 100.0);
        let camera_to_screen = self.camera.perspective_raw(aspect_ratio, 0.01, 100.0);

        let world_to_camera = self.camera.world_to_camera();

        let mut z_buffer = vec![vec![f64::INFINITY; width as usize + 1]; height as usize + 1];

        for shape in self.shapes.iter() {
            let shape_to_world = shape.shape_to_world();
            let overall = camera_to_screen * world_to_camera * shape_to_world;

            for tri in shape.triangles.iter() {
                let (a, b, c) = (
                    shape.vertices[tri[0]],
                    shape.vertices[tri[1]],
                    shape.vertices[tri[2]],
                );
                let a = overall.transform_point(&a);
                let b = overall.transform_point(&b);
                let c = overall.transform_point(&c);
                if a.x < -aspect_ratio || a.x > aspect_ratio || a.y < -1.0 || a.y > 1.0 {
                    continue;
                }
                if b.x < -aspect_ratio || b.x > aspect_ratio || b.y < -1.0 || b.y > 1.0 {
                    continue;
                }
                if c.x < -aspect_ratio || c.x > aspect_ratio || c.y < -1.0 || c.y > 1.0 {
                    continue;
                }

                let z = (a.z + b.z + c.z) / 3.0;
                if z > 1.0 {
                    continue;
                }

                // println!("{}", t.z);

                // Ensure counterclockwise winding
                let (a, b, c) = sort_points_ccw((a, b, c));

                let ax = ((a.x + 1.0) * 0.5 * fwidth) as i32;
                let ay = ((-a.y + 1.0) * 0.5 * fwidth) as i32; // in calculation +Y up, but on screen +Y down
                let ax = ax.clamp(0, width as i32);
                let ay = ay.clamp(0, height as i32);
                let a = Pt3::new(ax as f64, ay as f64, a.z);
                let a_screen = embedded_graphics::geometry::Point::new(ax, ay);

                let bx = ((b.x + 1.0) * 0.5 * fwidth) as i32;
                let by = ((-b.y + 1.0) * 0.5 * fwidth) as i32; // in calculation +Y up, but on screen +Y down
                let bx = bx.clamp(0, width as i32);
                let by = by.clamp(0, height as i32);
                let b = Pt3::new(bx as f64, by as f64, b.z);
                let b_screen = embedded_graphics::geometry::Point::new(bx, by);

                let cx = ((c.x + 1.0) * 0.5 * fwidth) as i32;
                let cy = ((-c.y + 1.0) * 0.5 * fwidth) as i32; // in calculation +Y up, but on screen +Y down
                let cx = cx.clamp(0, width as i32);
                let cy = cy.clamp(0, height as i32);
                let c = Pt3::new(cx as f64, cy as f64, c.z);
                let c_screen = embedded_graphics::geometry::Point::new(cx, cy);

                let sxmin = a_screen.x.min(b_screen.x).min(c_screen.x);
                let symin = a_screen.y.min(b_screen.y).min(c_screen.y);
                let sxmax = a_screen.x.max(b_screen.x).max(c_screen.x);
                let symax = a_screen.y.max(b_screen.y).max(c_screen.y);

                match shape.style {
                    RenderingStyle::Wireframe { line_style } => {
                        embedded_graphics::primitives::Triangle::new(a_screen, b_screen, c_screen)
                            .into_styled(line_style)
                            .draw(frame)?
                    }

                    RenderingStyle::UnlitRandom => {
                        let color = Rgb888::new(rng.gen(), rng.gen(), rng.gen());
                        let tri = Triangle::new(a_screen, b_screen, c_screen);
                        //.into_styled(PrimitiveStyleBuilder::new().fill_color(color).build())
                        //.draw(frame)?;

                        let triangle_area = edge_function(a, b, c);

                        for pt_screen in tri.points() {
                            let Point { x, y } = pt_screen;
                            let pt = Pt3::new(x as f64, y as f64, 0.0);
                            // Barycentric coordinates to check if point inside triangle
                            let w0 = edge_function(b, c, pt) / triangle_area;
                            let w1 = edge_function(c, a, pt) / triangle_area;
                            let w2 = edge_function(a, b, pt) / triangle_area;

                            if is_pixel_overlapping((a, b, c), (w0, w1, w2)) {
                                // Interpolate Z coordinate
                                let z = 1.0 / (1.0 / a.z * w0 + 1.0 / b.z * w1 + 1.0 / c.z * w2);

                                // if z > 1.0 {
                                //     continue;
                                // }
                                if z < z_buffer[y as usize][x as usize] {
                                    z_buffer[y as usize][x as usize] = z;
                                    let mut p = embedded_graphics::prelude::Pixel::default();
                                    p.0 = pt_screen;
                                    p.1 = color;
                                    p.draw(frame)?;
                                }
                            }
                        }
                    }
                }
            }
        }

        // Draw the coordinates of the camera in the top left corner.
        let text = format!(
            "{:.2} {:.2} {:.2}",
            self.camera.location.x, self.camera.location.y, self.camera.location.z
        );
        let text = embedded_graphics::text::Text::new(
            &text,
            Point::new(0, 15),
            MonoTextStyle::new(&FONT_10X20, Rgb888::new(255, 255, 255)),
        );
        text.draw(frame)?;

        // Draw the rotation.
        let angles = self.camera.rotation.euler_angles();
        let text = format!(
            "R:{:.2} P:{:.2} Y:{:.2}",
            f64::to_degrees(angles.0),
            f64::to_degrees(angles.1),
            f64::to_degrees(angles.2)
        );
        let text = embedded_graphics::text::Text::new(
            &text,
            Point::new(0, 240),
            MonoTextStyle::new(&FONT_10X20, Rgb888::new(255, 255, 255)),
        );
        text.draw(frame)?;

        Ok(())
    }

    pub fn update(&mut self) {
        // self.camera.location += Vec3::new(0.01, -0.01, 0.01);
        self.shapes[0].rotation =
            self.shapes[0].rotation * Quat::from_axis_angle(&Vec3::y_axis(), 0.01);
    }
}
