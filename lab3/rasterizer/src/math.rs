use std::f64::consts::TAU;

use nalgebra::{Matrix2, Matrix4, UnitQuaternion, Vector3};

pub type Vec3 = Vector3<f64>;
pub type Quat = UnitQuaternion<f64>;
pub type Mat4 = Matrix4<f64>;

#[rustfmt::skip]
pub fn xy_rotation_mat(angle_rad: f64) -> Mat4 {
    let (s, c) = angle_rad.sin_cos();
    let o = 0.0;
    let i = 1.0;
    Mat4::new(
        c, -s, o, o,
        s, c,  o, o,
        o, o,  i, o,
        o, o,  o, i
    )
}

#[rustfmt::skip]
pub fn yz_rotation_mat(angle_rad: f64) -> Mat4 {
    let (s, c) = angle_rad.sin_cos();
    let o = 0.0;
    let i = 1.0;
    Mat4::new(
        i, o, o,  o,
        o, c, -s, o,
        o, s, c,  o,
        o, o, o,  i
    )
}

#[rustfmt::skip]
pub fn xz_rotation_mat(angle_rad: f64) -> Mat4 {
    let (s, c) = angle_rad.sin_cos();
    let o = 0.0;
    let i = 1.0;
    Mat4::new(
        c, o, -s, o,
        o, i, o,  o,
        s, o, c,  o,
        o, o, o,  i
    )
}

pub fn quat_to_rotmat(quat: Quat) -> Mat4 {
    let (xz, yz, xy) = quat.euler_angles();
    xz_rotation_mat(xz) * yz_rotation_mat(yz) * xy_rotation_mat(xy)
}

pub fn sort_points_clockwise(src: (Vec3, Vec3, Vec3)) -> (Vec3, Vec3, Vec3) {
    let (a, b, c) = src;
    let centroid = (a + b + c) / 3.0;
    let mut verts = [a, b, c];
    verts.sort_by(|a, b| {
        (f64::atan2(a[0] - centroid[0], a[1] - centroid[1]) % TAU)
            .partial_cmp(&(f64::atan2(b[0] - centroid[0], b[1] - centroid[1]) % TAU))
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    (verts[0], verts[1], verts[2])
}

pub fn edge_function(v0: Vec3, v1: Vec3, point: Vec3) -> f64 {
    Matrix2::new(
        point[0] - v0[0],
        point[1] - v0[1],
        v1[0] - v0[0],
        v1[1] - v0[1],
    )
    .determinant()
}

pub fn is_pixel_overlapping(
    screen_verts: (Vec3, Vec3, Vec3),
    barycentric: (f64, f64, f64),
) -> bool {
    let (a, b, c) = screen_verts;
    let edge0 = c - b;
    let edge1 = b - a;
    let edge2 = a - b;

    let o = 0.0;

    if barycentric.0 == o {
        if !(edge0[1] == o && edge0[0] > o || edge0[1] > o) {
            return false;
        }
    } else if barycentric.0 <= o {
        return false;
    }

    if barycentric.1 == o {
        if !(edge1[1] == o && edge1[0] > o || edge1[1] > o) {
            return false;
        }
    } else if barycentric.1 <= o {
        return false;
    }

    if barycentric.2 == o {
        if !(edge2[1] == o && edge2[0] > o || edge2[1] > o) {
            return false;
        }
    } else if barycentric.2 <= o {
        return false;
    }

    true
}
