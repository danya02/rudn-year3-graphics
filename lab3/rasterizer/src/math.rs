use std::f64::consts::TAU;

use nalgebra::{Matrix2, Matrix4, Point3, UnitQuaternion, Vector3};

pub type Vec3 = Vector3<f64>;
pub type Pt3 = Point3<f64>;
pub type Quat = UnitQuaternion<f64>;
pub type Mat4 = Matrix4<f64>;

pub fn sort_points_ccw(src: (Pt3, Pt3, Pt3)) -> (Pt3, Pt3, Pt3) {
    let (a, b, c) = src;
    let centroid = (a.coords + b.coords + c.coords) / 3.0;
    let mut verts = [a, b, c];
    verts.sort_by(|a, b| {
        (f64::atan2(a[0] - centroid[0], a[1] - centroid[1]) % TAU)
            .partial_cmp(&(f64::atan2(b[0] - centroid[0], b[1] - centroid[1]) % TAU))
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    (verts[0], verts[1], verts[2])
}

pub fn edge_function(v0: Pt3, v1: Pt3, point: Pt3) -> f64 {
    Matrix2::new(point.x - v0.x, point.y - v0.y, v1.x - v0.x, v1.y - v0.y).determinant()
}

pub fn is_pixel_overlapping(screen_verts: (Pt3, Pt3, Pt3), barycentric: (f64, f64, f64)) -> bool {
    let (a, b, c) = screen_verts;
    let edge0 = c - b;
    let edge1 = a - b;
    let edge2 = b - a;

    let o = 0.0;

    if if barycentric.0 == o {
        edge0.y == o && edge0.x > o || edge0.y > o
    } else {
        barycentric.0 <= o
    } {
        return false;
    }

    if if barycentric.1 == o {
        edge1.y == o && edge1.x > o || edge1.y > o
    } else {
        barycentric.1 <= o
    } {
        return false;
    }

    if if barycentric.2 == o {
        edge2.y == o && edge2.x > o || edge2.y > o
    } else {
        barycentric.2 <= o
    } {
        return false;
    }

    true
}
