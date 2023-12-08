use std::{collections::HashMap, io::Cursor};

use embedded_graphics::pixelcolor::Rgb888;

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
        style: crate::world::RenderingStyle::GlobalLight {
            color: Rgb888::new(0, 0, 255),
            line_color: Some(Rgb888::new(0, 0, 255)),
        },
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
            .map(|v| Pt3::new(v[0] as f64, v[2] as f64, v[1] as f64))
            .collect(),
        triangles: stl.faces.iter().map(|f| f.vertices).collect(),
        style: crate::world::RenderingStyle::GlobalLight {
            color: Rgb888::new(255, 0, 0),
            line_color: None,
        },
    }
}

const CUBE_DATA: &'static [u8] = include_bytes!("cube.stl");

pub fn get_cube() -> Shape {
    let stl = stl_io::read_stl(&mut Cursor::new(CUBE_DATA)).unwrap();
    Shape {
        origin: Pt3::origin(),
        rotation: Quat::identity(),
        scale_factor: 1.0,
        vertices: stl
            .vertices
            .iter()
            .map(|v| Pt3::new(v[0] as f64, v[2] as f64, v[1] as f64))
            .collect(),
        triangles: stl.faces.iter().map(|f| f.vertices).collect(),
        style: crate::world::RenderingStyle::GlobalLight {
            color: Rgb888::new(255, 0, 255),
            line_color: None,
        },
    }
}

const TETRAHEDRON_DATA: &'static [u8] = include_bytes!("tetrahedron.stl");

pub fn get_tetrahedron() -> Shape {
    let stl = stl_io::read_stl(&mut Cursor::new(TETRAHEDRON_DATA)).unwrap();
    Shape {
        origin: Pt3::origin(),
        rotation: Quat::identity(),
        scale_factor: 1.0,
        vertices: stl
            .vertices
            .iter()
            .map(|v| Pt3::new(v[0] as f64, v[2] as f64, v[1] as f64))
            .collect(),
        triangles: stl.faces.iter().map(|f| f.vertices).collect(),
        style: crate::world::RenderingStyle::GlobalLight {
            color: Rgb888::new(0, 255, 0),
            line_color: None,
        },
    }
}

const ICOSAHEDRON_DATA: &'static [u8] = include_bytes!("icosahedron.stl");

pub fn get_icosahedron() -> Shape {
    let stl = stl_io::read_stl(&mut Cursor::new(ICOSAHEDRON_DATA)).unwrap();
    Shape {
        origin: Pt3::origin(),
        rotation: Quat::identity(),
        scale_factor: 1.0,
        vertices: stl
            .vertices
            .iter()
            .map(|v| Pt3::new(v[0] as f64, v[2] as f64, v[1] as f64))
            .collect(),
        triangles: stl.faces.iter().map(|f| f.vertices).collect(),
        style: crate::world::RenderingStyle::GlobalLight {
            color: Rgb888::new(255, 255, 0),
            line_color: None,
        },
    }
}

pub fn get_floor() -> Vec<Shape> {
    let mut output = vec![];
    // For each of the squares, output a shape with two triangles.
    for x in -10..10 {
        for y in -10..10 {
            let shape = Shape {
                origin: Pt3::new(x as f64, 0.0, y as f64),
                rotation: Quat::identity(),
                scale_factor: 1.0,
                vertices: vec![
                    Pt3::new(0.0, 0.0, 0.0),
                    Pt3::new(1.0, 0.0, 0.0),
                    Pt3::new(0.0, 0.0, 1.0),
                    Pt3::new(1.0, 0.0, 1.0),
                ],
                triangles: vec![[0, 1, 2], [2, 3, 1]],
                style: crate::world::RenderingStyle::GlobalLight {
                    color: Rgb888::new(64, 0, 0),
                    line_color: None,
                },
            };

            output.push(shape);
        }
    }
    output
}

pub fn get_axes() -> Vec<Shape> {
    let mut output = vec![];
    // X
    let shape = Shape {
        origin: Pt3::new(0.0, 0.1, 0.0),
        rotation: Quat::identity(),
        scale_factor: 1.0,
        vertices: vec![
            Pt3::new(0.0, 0.0, 0.0),
            Pt3::new(10.0, 0.0, 0.0),
            Pt3::new(0.0, 0.0, 0.1),
            Pt3::new(10.0, 0.0, 0.1),
        ],
        triangles: vec![[0, 1, 2], [2, 3, 1]],
        style: crate::world::RenderingStyle::Unlit {
            color: Rgb888::new(255, 0, 0),
            line_color: None,
        },
    };

    output.push(shape);

    // Y
    let shape = Shape {
        origin: Pt3::new(0.0, 0.1, 0.0),
        rotation: Quat::identity(),
        scale_factor: 1.0,
        vertices: vec![
            Pt3::new(0.0, 0.0, 0.0),
            Pt3::new(0.0, 10.0, 0.0),
            Pt3::new(0.0, 0.0, 0.1),
            Pt3::new(0.0, 10.0, 0.1),
        ],
        triangles: vec![[0, 1, 2], [2, 3, 1]],
        style: crate::world::RenderingStyle::Unlit {
            color: Rgb888::new(0, 255, 0),
            line_color: None,
        },
    };

    output.push(shape);

    // Z
    let shape = Shape {
        origin: Pt3::new(0.0, 0.1, 0.0),
        rotation: Quat::identity(),
        scale_factor: 1.0,
        vertices: vec![
            Pt3::new(0.0, 0.0, 0.0),
            Pt3::new(0.0, 0.1, 0.0),
            Pt3::new(0.0, 0.0, 10.0),
            Pt3::new(0.0, 0.1, 10.0),
        ],
        triangles: vec![[0, 1, 2], [2, 3, 1]],
        style: crate::world::RenderingStyle::Unlit {
            color: Rgb888::new(0, 0, 255),
            line_color: None,
        },
    };

    output.push(shape);

    output
}
