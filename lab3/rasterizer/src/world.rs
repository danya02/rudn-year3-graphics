use std::convert::Infallible;

use bresenham::Bresenham;
use embedded_graphics::{
    draw_target::DrawTarget,
    geometry::OriginDimensions,
    pixelcolor::Rgb888,
    primitives::{Primitive, PrimitiveStyleBuilder},
    Drawable,
};
use nalgebra::{
    Isometry3, Matrix, Matrix4, OPoint, Orthographic3, Perspective3, Point3, Rotation3, Scale3,
    Translation3, UnitQuaternion, Vector2, Vector3, Vector4,
};

use crate::math::{
    edge_function, is_pixel_overlapping, quat_to_rotmat, sort_points_clockwise,
    transform_point_by_projection, Mat4, Pt3, Quat, Vec3,
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
    pub triangles: Vec<(usize, usize, usize)>,
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
}

impl World {
    pub fn rasterize<D: DrawTarget<Color = Rgb888, Error = Infallible> + OriginDimensions>(
        &self,
        frame: &mut D,
    ) -> Result<(), Infallible> {
        frame.clear(Rgb888::new(0, 0, 0))?;

        let (mut maxx, mut maxy, mut maxz) =
            (f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        let (mut minx, mut miny, mut minz) = (f64::INFINITY, f64::INFINITY, f64::INFINITY);
        for shape in self.shapes.iter() {
            let shape_to_world = shape.shape_to_world();
            for tri in shape.triangles.iter() {
                for point in [
                    shape.vertices[tri.0],
                    shape.vertices[tri.1],
                    shape.vertices[tri.2],
                ]
                .iter()
                {
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
        let min = Pt3::new(minx, miny, minz);
        let max = Pt3::new(maxx, maxy, maxz);
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

        let aspect_ratio = fwidth / fheight;
        let r = max * aspect_ratio;
        let l = -r;
        let t = max;
        let b = -t;

        // let camera_to_screen = self.camera.ortho_raw(r, l, t, b, 0.01, 100.0);
        let camera_to_screen = self.camera.perspective_raw(aspect_ratio, 0.01, 100.0);

        let world_to_camera = self.camera.world_to_camera();

        for shape in self.shapes.iter() {
            let shape_to_world = shape.shape_to_world();

            for tri in shape.triangles.iter() {
                for point in [
                    shape.vertices[tri.0],
                    shape.vertices[tri.1],
                    shape.vertices[tri.2],
                ]
                .iter()
                {
                    let t = point;
                    let t = shape_to_world.transform_point(&t);
                    let t = world_to_camera.transform_point(&t);
                    let t = camera_to_screen.transform_point(&t);
                    // if t.x < -aspect_ratio || t.x > aspect_ratio || t.y < -1.0 || t.y > 1.0 {
                    //     continue;
                    // }
                    // println!("{}", t.z);

                    let tx = ((t.x + 1.0) * 0.5 * fwidth) as i32;
                    let ty = ((-t.y + 1.0) * 0.5 * fwidth) as i32; // in calculation +Y up, but on screen +Y down
                    let tx = tx.clamp(0, width as i32);
                    let ty = ty.clamp(0, height as i32);
                    let t_screen = embedded_graphics::geometry::Point::new(tx, ty);
                    embedded_graphics::primitives::Rectangle::new(
                        t_screen,
                        embedded_graphics::geometry::Size::new(3, 3),
                    )
                    .into_styled(
                        PrimitiveStyleBuilder::new()
                            .fill_color(if t.z < 0.0 {
                                Rgb888::new(0, 0, 255)
                            } else {
                                Rgb888::new(255, 0, 0)
                            })
                            .build(),
                    )
                    .draw(frame)?;
                }
            }
        }
        Ok(())
    }

    pub fn update(&mut self) {
        // self.camera.location += Vec3::new(0.01, -0.01, 0.01);
        self.shapes[0].rotation =
            self.shapes[0].rotation * Quat::from_axis_angle(&Vec3::y_axis(), 0.01);
    }
}
