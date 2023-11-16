use std::collections::HashMap;

use crate::{
    math::{Quat, Vec3},
    world::Shape,
};

const TEAPOT_TRIS: &'static str = include_str!("./teapot.tris");

pub fn get_teapot() -> Shape {
    let mut verts = vec![];
    let mut vert_idxs = HashMap::new();
    let mut tris = vec![];
    let mut words = TEAPOT_TRIS.split_whitespace();

    let tri_count = words.next().unwrap().parse().unwrap();
    for _ in 0..tri_count {
        let a = {
            let x = words.next().unwrap();
            let y = words.next().unwrap();
            let z = words.next().unwrap();
            if let Some(idx) = vert_idxs.get(&(x, y, z)) {
                *idx
            } else {
                let i = verts.len();
                vert_idxs.insert((x, y, z), i);

                let x = x.parse().unwrap();
                let y = y.parse().unwrap();
                let z = z.parse().unwrap();
                verts.push(Vec3::new(x, y, z));
                i
            }
        };
        let b = {
            let x = words.next().unwrap();
            let y = words.next().unwrap();
            let z = words.next().unwrap();
            if let Some(idx) = vert_idxs.get(&(x, y, z)) {
                *idx
            } else {
                let i = verts.len();
                vert_idxs.insert((x, y, z), i);

                let x = x.parse().unwrap();
                let y = y.parse().unwrap();
                let z = z.parse().unwrap();
                verts.push(Vec3::new(x, y, z));
                i
            }
        };
        let c = {
            let x = words.next().unwrap();
            let y = words.next().unwrap();
            let z = words.next().unwrap();

            if let Some(idx) = vert_idxs.get(&(x, y, z)) {
                *idx
            } else {
                let i = verts.len();
                vert_idxs.insert((x, y, z), i);

                let x = x.parse().unwrap();
                let y = y.parse().unwrap();
                let z = z.parse().unwrap();
                verts.push(Vec3::new(x, y, z));
                i
            }
        };

        tris.push((a, b, c));
    }

    Shape {
        origin: Vec3::zeros(),
        rotation: Quat::identity(),
        vertices: verts,
        triangles: tris,
    }
}
