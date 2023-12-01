use std::{collections::HashMap, io::Cursor};

use embedded_graphics::{pixelcolor::Rgb888, primitives::PrimitiveStyleBuilder};

use crate::{
    math::{Pt3, Quat},
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
                verts.push(Pt3::new(x, y, z));
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
                verts.push(Pt3::new(x, y, z));
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
                verts.push(Pt3::new(x, y, z));
                i
            }
        };

        tris.push([a, b, c]);
    }

    Shape {
        origin: Pt3::origin(),
        rotation: Quat::identity(),
        vertices: verts,
        triangles: tris,
        scale_factor: 1.0,
        style: crate::world::RenderingStyle::UnlitRandom,
        // style: crate::world::RenderingStyle::Wireframe {
        //     line_style: PrimitiveStyleBuilder::new()
        //         .stroke_width(2)
        //         .stroke_color(Rgb888::new(255, 0, 0))
        //         .build(),
        // },
    }
}

const SUZANNE_DATA: &'static [u8] = include_bytes!("suzanne.stl");

pub fn get_suzanne() -> Shape {
    let stl = stl_io::read_stl(&mut Cursor::new(SUZANNE_DATA)).unwrap();
    Shape {
        origin: Pt3::origin(),
        rotation: Quat::identity(),
        scale_factor: 1.0,
        vertices: stl
            .vertices
            .iter()
            .map(|v| Pt3::new(v[0] as f64, v[1] as f64, v[2] as f64))
            .collect(),
        triangles: stl.faces.iter().map(|f| f.vertices).collect(),
        style: crate::world::RenderingStyle::Wireframe {
            line_style: PrimitiveStyleBuilder::new()
                .stroke_width(2)
                .stroke_color(Rgb888::new(255, 0, 0))
                .build(),
        },
    }
}
