use std::{convert::Infallible, f64::consts::FRAC_PI_4};

use bresenham::Bresenham;
use embedded_graphics::{
    draw_target::DrawTarget,
    geometry::{OriginDimensions, Point},
    mono_font::{ascii::FONT_10X20, MonoTextStyle},
    pixelcolor::{Rgb888, RgbColor},
    primitives::{PointsIter, Primitive, PrimitiveStyle, Triangle},
    Drawable,
};
use nalgebra::{
    Matrix4, Orthographic3, Perspective3, Quaternion, Rotation3, Scale3, Translation3,
    UnitQuaternion, UnitVector3,
};
use rand::{Rng, SeedableRng};

use crate::{
    math::{edge_function, is_pixel_overlapping, sort_points_ccw, Mat4, Pt3, Quat, Vec3},
    PixelCounter,
};

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

    pub light_vector: UnitVector3<f64>,
    pub light_vector_eulers: Vec3,
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

    pub fn recalculate_light(&mut self) {
        let rot = UnitQuaternion::from_euler_angles(
            self.light_vector_eulers.x,
            self.light_vector_eulers.y,
            self.light_vector_eulers.z,
        );
        self.light_vector = UnitVector3::new_normalize(rot.transform_vector(&-Vec3::y_axis()));
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
    /// However, filling the triangle will result in faces being shown in an unnatural ordering, ignoring the Z-buffer.
    Wireframe { line_style: PrimitiveStyle<Rgb888> },

    /// Draw each face in a random color, and use the Z-buffer for rendering.
    UnlitRandom,

    /// Draw every face with the same color, and use the Z-buffer for rendering.
    Unlit {
        color: Rgb888,
        line_color: Option<Rgb888>,
    },

    /// Draw each face in a color, dimmed by the dot product of the view * normal vector.
    ViewLight {
        color: Rgb888,
        line_color: Option<Rgb888>,
    },

    /// Draw each face in a color, dimmed by the dot product of the global light * normal vector.
    GlobalLight {
        color: Rgb888,
        line_color: Option<Rgb888>,
    },
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

    pub fn with_scale(mut self, scale: f64) -> Self {
        self.scale_factor = scale;
        self
    }
}

impl World {
    pub fn rasterize<
        D: DrawTarget<Color = Rgb888, Error = Infallible> + OriginDimensions + PixelCounter,
    >(
        &mut self,
        frame: &mut D,
    ) -> Result<(), Infallible> {
        frame.clear(Rgb888::new(0, 0, 0))?;
        frame.reset_counter();

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

        let point_on_screen = |x, y| x > 0 && y > 0 && x < iwidth && y < iheight;

        let aspect_ratio = fwidth / fheight;
        let r = max * aspect_ratio;
        let l = -r;
        let t = max;
        let b = -t;

        let camera_to_screen = match self.camera.projection {
            Projection::Perspective => self.camera.perspective_raw(aspect_ratio, 0.01, 100.0),
            Projection::Orthographic => {
                let aspect_ratio = fwidth / fheight;
                let fac = self.camera.fov_radians / FRAC_PI_4;
                let r = max * aspect_ratio * fac;
                let l = -r;
                let t = max * fac;
                let b = -t;
                self.camera.ortho_raw(r, l, t, b, 0.01, 100.0)
            }
        };

        let world_to_camera = self.camera.world_to_camera();

        let mut z_buffer = vec![vec![1.0; width as usize + 1]; height as usize + 1];
        let mut drawn_triangles = 0;

        let view_vector = world_to_camera
            .try_inverse()
            .unwrap()
            .transform_vector(&UnitVector3::new_unchecked(Vec3::y()));
        let view_vector = UnitVector3::new_normalize(view_vector.normalize());

        for shape in self.shapes.iter() {
            let shape_to_world = shape.shape_to_world();
            let overall = camera_to_screen * world_to_camera * shape_to_world;

            for tri in shape.triangles.iter() {
                let (world_a, world_b, world_c) = (
                    shape.vertices[tri[0]],
                    shape.vertices[tri[1]],
                    shape.vertices[tri[2]],
                );
                let a = overall.transform_point(&world_a);
                let b = overall.transform_point(&world_b);
                let c = overall.transform_point(&world_c);

                let min_x = a.x.min(b.x).min(c.x);
                let max_x = a.x.max(b.x).max(c.x);
                let min_y = a.y.min(b.y).min(c.y);
                let max_y = a.y.max(b.y).max(c.y);

                let buf = 1.5;
                if min_x < -aspect_ratio * buf
                    || max_x > aspect_ratio * buf
                    || min_y < -1.0 * buf
                    || max_y > 1.0 * buf
                {
                    continue;
                }

                let z = a.z.min(b.z).min(c.z);
                if z > 1.0 {
                    continue;
                }

                // println!("{}", t.z);

                // Ensure counterclockwise winding
                let (a, b, c) = sort_points_ccw((a, b, c));

                let ax = ((a.x + 1.0) * 0.5 * fwidth).round() as i32;
                let ay = ((-a.y + 1.0) * 0.5 * fwidth).round() as i32; // in calculation +Y up, but on screen +Y down
                                                                       // let ax = ax.clamp(0, width as i32);
                                                                       // let ay = ay.clamp(0, height as i32);
                let a = Pt3::new(ax as f64, ay as f64, a.z);
                let a_screen = embedded_graphics::geometry::Point::new(ax, ay);

                let bx = ((b.x + 1.0) * 0.5 * fwidth).round() as i32;
                let by = ((-b.y + 1.0) * 0.5 * fwidth).round() as i32; // in calculation +Y up, but on screen +Y down
                                                                       // let bx = bx.clamp(0, width as i32);
                                                                       // let by = by.clamp(0, height as i32);
                let b = Pt3::new(bx as f64, by as f64, b.z);
                let b_screen = embedded_graphics::geometry::Point::new(bx, by);

                let cx = ((c.x + 1.0) * 0.5 * fwidth).round() as i32;
                let cy = ((-c.y + 1.0) * 0.5 * fwidth).round() as i32; // in calculation +Y up, but on screen +Y down
                                                                       // let cx = cx.clamp(0, width as i32);
                                                                       // let cy = cy.clamp(0, height as i32);
                let c = Pt3::new(cx as f64, cy as f64, c.z);
                let c_screen = embedded_graphics::geometry::Point::new(cx, cy);

                // let sxmin = a_screen.x.min(b_screen.x).min(c_screen.x);
                // let symin = a_screen.y.min(b_screen.y).min(c_screen.y);
                // let sxmax = a_screen.x.max(b_screen.x).max(c_screen.x);
                // let symax = a_screen.y.max(b_screen.y).max(c_screen.y);

                drawn_triangles += 1;

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
                            if !point_on_screen(x, y) {
                                continue;
                            }

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

                    RenderingStyle::Unlit { color, line_color } => {
                        let tri = Triangle::new(a_screen, b_screen, c_screen);
                        //.into_styled(PrimitiveStyleBuilder::new().fill_color(color).build())
                        //.draw(frame)?;

                        let triangle_area = edge_function(a, b, c);

                        for pt_screen in tri.points() {
                            let Point { x, y } = pt_screen;
                            if !point_on_screen(x, y) {
                                continue;
                            }

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

                        if let Some(color) = line_color {
                            for (x, y) in Bresenham::new(
                                (tri.vertices[0].x as isize, tri.vertices[0].y as isize),
                                (tri.vertices[1].x as isize, tri.vertices[1].y as isize),
                            )
                            .chain(Bresenham::new(
                                (tri.vertices[0].x as isize, tri.vertices[0].y as isize),
                                (tri.vertices[2].x as isize, tri.vertices[2].y as isize),
                            ))
                            .chain(Bresenham::new(
                                (tri.vertices[1].x as isize, tri.vertices[1].y as isize),
                                (tri.vertices[2].x as isize, tri.vertices[2].y as isize),
                            )) {
                                if !point_on_screen(x as i32, y as i32) {
                                    continue;
                                }
                                let pt = Pt3::new(x as f64, y as f64, 0.0);
                                // Barycentric coordinates to check if point inside triangle
                                let w0 = edge_function(b, c, pt) / triangle_area;
                                let w1 = edge_function(c, a, pt) / triangle_area;
                                let w2 = edge_function(a, b, pt) / triangle_area;
                                // if is_pixel_overlapping((a, b, c), (w0, w1, w2)) {
                                // Interpolate Z coordinate
                                let z = 1.0 / (1.0 / a.z * w0 + 1.0 / b.z * w1 + 1.0 / c.z * w2);

                                // if z > 1.0 {
                                //     continue;
                                // }
                                if z < z_buffer[y as usize][x as usize] {
                                    z_buffer[y as usize][x as usize] = z;
                                    let mut p = embedded_graphics::prelude::Pixel::default();
                                    p.0 = Point::new(x as i32, y as i32);
                                    p.1 = color;
                                    p.draw(frame)?;
                                }
                                //}
                            }
                        }
                    }

                    RenderingStyle::ViewLight { color, line_color } => {
                        let rcolor = color.r();
                        let gcolor = color.g();
                        let bcolor = color.b();
                        let tri = Triangle::new(a_screen, b_screen, c_screen);
                        //.into_styled(PrimitiveStyleBuilder::new().fill_color(color).build())
                        //.draw(frame)?;

                        let normal = UnitVector3::new_normalize(
                            (world_c - world_a).cross(&(world_b - world_a)),
                        );

                        let dot = normal.dot(&view_vector) as f32;
                        let dot = dot.abs();
                        let dot = dot.clamp(0.0, 1.0);
                        // if dot.is_nan() {
                        //     log::warn!("{world_a:?} {world_b:?}");
                        // }

                        let color = Rgb888::new(
                            (rcolor as f32 * dot) as u8,
                            (gcolor as f32 * dot) as u8,
                            (bcolor as f32 * dot) as u8,
                        );

                        let triangle_area = edge_function(a, b, c);

                        for pt_screen in tri.points() {
                            let Point { x, y } = pt_screen;
                            if !point_on_screen(x, y) {
                                continue;
                            }

                            let pt = Pt3::new(x as f64, y as f64, 0.0);
                            // Barycentric coordinates to check if point inside triangle
                            let w0 = edge_function(b, c, pt) / triangle_area;
                            let w1 = edge_function(c, a, pt) / triangle_area;
                            let w2 = edge_function(a, b, pt) / triangle_area;

                            //if is_pixel_overlapping((a, b, c), (w0, w1, w2)) {
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
                            //}
                        }

                        if let Some(color) = line_color {
                            for (x, y) in Bresenham::new(
                                (tri.vertices[0].x as isize, tri.vertices[0].y as isize),
                                (tri.vertices[1].x as isize, tri.vertices[1].y as isize),
                            )
                            .chain(Bresenham::new(
                                (tri.vertices[0].x as isize, tri.vertices[0].y as isize),
                                (tri.vertices[2].x as isize, tri.vertices[2].y as isize),
                            ))
                            .chain(Bresenham::new(
                                (tri.vertices[1].x as isize, tri.vertices[1].y as isize),
                                (tri.vertices[2].x as isize, tri.vertices[2].y as isize),
                            )) {
                                if !point_on_screen(x as i32, y as i32) {
                                    continue;
                                }
                                let pt = Pt3::new(x as f64, y as f64, 0.0);
                                // Barycentric coordinates to check if point inside triangle
                                let w0 = edge_function(b, c, pt) / triangle_area;
                                let w1 = edge_function(c, a, pt) / triangle_area;
                                let w2 = edge_function(a, b, pt) / triangle_area;
                                // if is_pixel_overlapping((a, b, c), (w0, w1, w2)) {
                                // Interpolate Z coordinate
                                let z = 1.0 / (1.0 / a.z * w0 + 1.0 / b.z * w1 + 1.0 / c.z * w2);

                                // if z > 1.0 {
                                //     continue;
                                // }
                                if z < z_buffer[y as usize][x as usize] {
                                    z_buffer[y as usize][x as usize] = z;
                                    let mut p = embedded_graphics::prelude::Pixel::default();
                                    p.0 = Point::new(x as i32, y as i32);
                                    p.1 = color;
                                    p.draw(frame)?;
                                }
                                //}
                            }
                        }
                    }
                    RenderingStyle::GlobalLight { color, line_color } => {
                        let rcolor = color.r();
                        let gcolor = color.g();
                        let bcolor = color.b();
                        let tri = Triangle::new(a_screen, b_screen, c_screen);
                        //.into_styled(PrimitiveStyleBuilder::new().fill_color(color).build())
                        //.draw(frame)?;

                        let normal = UnitVector3::new_normalize(
                            (world_c - world_a).cross(&(world_b - world_a)),
                        );

                        let dot = normal.dot(&self.camera.light_vector) as f32;
                        let dot = dot.abs();
                        let dot = dot.clamp(0.0, 1.0);

                        let color = Rgb888::new(
                            (rcolor as f32 * dot) as u8,
                            (gcolor as f32 * dot) as u8,
                            (bcolor as f32 * dot) as u8,
                        );

                        let triangle_area = edge_function(a, b, c);

                        for pt_screen in tri.points() {
                            let Point { x, y } = pt_screen;
                            if !point_on_screen(x, y) {
                                continue;
                            }

                            let pt = Pt3::new(x as f64, y as f64, 0.0);
                            // Barycentric coordinates to check if point inside triangle
                            let w0 = edge_function(b, c, pt) / triangle_area;
                            let w1 = edge_function(c, a, pt) / triangle_area;
                            let w2 = edge_function(a, b, pt) / triangle_area;

                            //if is_pixel_overlapping((a, b, c), (w0, w1, w2)) {
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
                            //}
                        }

                        if let Some(color) = line_color {
                            for (x, y) in Bresenham::new(
                                (tri.vertices[0].x as isize, tri.vertices[0].y as isize),
                                (tri.vertices[1].x as isize, tri.vertices[1].y as isize),
                            )
                            .chain(Bresenham::new(
                                (tri.vertices[0].x as isize, tri.vertices[0].y as isize),
                                (tri.vertices[2].x as isize, tri.vertices[2].y as isize),
                            ))
                            .chain(Bresenham::new(
                                (tri.vertices[1].x as isize, tri.vertices[1].y as isize),
                                (tri.vertices[2].x as isize, tri.vertices[2].y as isize),
                            )) {
                                if !point_on_screen(x as i32, y as i32) {
                                    continue;
                                }
                                let pt = Pt3::new(x as f64, y as f64, 0.0);
                                // Barycentric coordinates to check if point inside triangle
                                let w0 = edge_function(b, c, pt) / triangle_area;
                                let w1 = edge_function(c, a, pt) / triangle_area;
                                let w2 = edge_function(a, b, pt) / triangle_area;
                                // if is_pixel_overlapping((a, b, c), (w0, w1, w2)) {
                                // Interpolate Z coordinate
                                let z = 1.0 / (1.0 / a.z * w0 + 1.0 / b.z * w1 + 1.0 / c.z * w2);

                                // if z > 1.0 {
                                //     continue;
                                // }
                                if z < z_buffer[y as usize][x as usize] {
                                    z_buffer[y as usize][x as usize] = z;
                                    let mut p = embedded_graphics::prelude::Pixel::default();
                                    p.0 = Point::new(x as i32, y as i32);
                                    p.1 = color;
                                    p.draw(frame)?;
                                }
                                //}
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
            Point::new(0, 480),
            MonoTextStyle::new(&FONT_10X20, Rgb888::new(255, 255, 255)),
        );
        text.draw(frame)?;
        let text = format!("Tris: {drawn_triangles}; pixels: {}", frame.get_counter());
        let text = embedded_graphics::text::Text::new(
            &text,
            Point::new(0, 460),
            MonoTextStyle::new(&FONT_10X20, Rgb888::new(255, 255, 255)),
        );
        text.draw(frame)?;

        let text = format!("FOV: {:.2}", self.camera.fov_radians.to_degrees());
        let text = embedded_graphics::text::Text::new(
            &text,
            Point::new(0, 440),
            MonoTextStyle::new(&FONT_10X20, Rgb888::new(255, 255, 255)),
        );
        text.draw(frame)?;

        let text = format!(
            "Global light: {:.2} {:.2} {:.2}",
            self.camera.light_vector.x, self.camera.light_vector.y, self.camera.light_vector.z
        );
        let text = embedded_graphics::text::Text::new(
            &text,
            Point::new(0, 420),
            MonoTextStyle::new(&FONT_10X20, Rgb888::new(255, 255, 255)),
        );
        text.draw(frame)?;

        Ok(())
    }

    pub fn update(&mut self) {
        // self.camera.location += Vec3::new(0.01, -0.01, 0.01);
        self.shapes[0].rotation =
            self.shapes[0].rotation * Quat::from_axis_angle(&Vec3::y_axis(), 0.01);

        self.shapes[1].rotation = Quat::look_at_rh(
            &(self.shapes[1].origin - self.camera.location).coords,
            &Vec3::y(),
        );

        self.shapes.last_mut().unwrap().origin = self.camera.location.into();
    }
}
