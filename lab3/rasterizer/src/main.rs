#![deny(clippy::all)]
#![forbid(unsafe_code)]

use embedded_graphics::draw_target::DrawTarget;
use embedded_graphics::geometry::OriginDimensions;
use embedded_graphics::pixelcolor::raw::ToBytes;
use embedded_graphics::primitives::PrimitiveStyleBuilder;
use error_iter::ErrorIter as _;
use log::error;
use math::{Pt3, Quat, Vec3};
use nalgebra::UnitQuaternion;
use pixels::{Pixels, SurfaceTexture};
use shapes::{
    get_axes, get_cube, get_floor, get_icosahedron, get_suzanne, get_teapot, get_tetrahedron,
};
use std::convert::Infallible;
use std::f64::consts::PI;
use std::rc::Rc;
use winit::dpi::{LogicalSize, Pixel};
use winit::event::{Event, VirtualKeyCode};
use winit::event_loop::{ControlFlow, EventLoop};
use winit::window::WindowBuilder;
use winit_input_helper::WinitInputHelper;
use world::Camera;

mod math;
mod shapes;
mod world;

const WIDTH: u32 = 640;
const HEIGHT: u32 = 480;

struct ArrayDrawTarget<'a> {
    pub array: &'a mut [u8],
    pub width: usize,
    pub pixels_drawn: usize,
}

impl<'a> OriginDimensions for ArrayDrawTarget<'a> {
    fn size(&self) -> embedded_graphics::prelude::Size {
        embedded_graphics::prelude::Size::new(
            self.width as u32,
            ((self.array.len() / 4) / self.width) as u32,
        )
    }
}

pub trait PixelCounter {
    fn get_counter(&self) -> usize;
    fn reset_counter(&mut self);
}

impl PixelCounter for ArrayDrawTarget<'_> {
    fn get_counter(&self) -> usize {
        self.pixels_drawn
    }

    fn reset_counter(&mut self) {
        self.pixels_drawn = 0;
    }
}

impl<'a> DrawTarget for ArrayDrawTarget<'a> {
    type Color = embedded_graphics::pixelcolor::Rgb888;

    type Error = Infallible;

    fn draw_iter<I>(&mut self, pixels: I) -> Result<(), Self::Error>
    where
        I: IntoIterator<Item = embedded_graphics::prelude::Pixel<Self::Color>>,
    {
        let size = self.size();
        for pixel in pixels {
            let pt = pixel.0;
            if pt.x >= size.width as i32 || pt.y >= size.height as i32 || pt.x < 0 || pt.y < 0 {
                continue;
            }
            let x = pt.x as usize;
            let y = pt.y as usize;

            let idx = 4 * (x + (self.width * y));
            // if idx + 3 > self.array.len() {
            //     continue;
            // }

            let color = &pixel.1.to_be_bytes();
            self.array[idx..idx + 3].copy_from_slice(color);
            self.pixels_drawn += 1;
        }

        Ok(())
    }
}

fn main() {
    #[cfg(target_arch = "wasm32")]
    {
        std::panic::set_hook(Box::new(console_error_panic_hook::hook));
        console_log::init_with_level(log::Level::Trace).expect("error initializing logger");

        wasm_bindgen_futures::spawn_local(run());
    }

    #[cfg(not(target_arch = "wasm32"))]
    {
        env_logger::init();

        pollster::block_on(run());
    }
}

async fn run() {
    let event_loop = EventLoop::new();
    let window = {
        let size = LogicalSize::new(WIDTH as f64, HEIGHT as f64);
        WindowBuilder::new()
            .with_title("Hello Pixels + Web")
            .with_inner_size(size)
            .with_min_inner_size(size)
            .build(&event_loop)
            .expect("WindowBuilder error")
    };

    let window = Rc::new(window);

    #[cfg(target_arch = "wasm32")]
    {
        use wasm_bindgen::JsCast;
        use winit::platform::web::WindowExtWebSys;

        // Retrieve current width and height dimensions of browser client window
        let get_window_size = || {
            let client_window = web_sys::window().unwrap();
            LogicalSize::new(
                client_window.inner_width().unwrap().as_f64().unwrap(),
                client_window.inner_height().unwrap().as_f64().unwrap(),
            )
        };

        let window = Rc::clone(&window);

        // Initialize winit window with current dimensions of browser client
        window.set_inner_size(get_window_size());

        let client_window = web_sys::window().unwrap();

        // Attach winit canvas to body element
        web_sys::window()
            .and_then(|win| win.document())
            .and_then(|doc| doc.body())
            .and_then(|body| {
                body.append_child(&web_sys::Element::from(window.canvas()))
                    .ok()
            })
            .expect("couldn't append canvas to document body");

        // Listen for resize event on browser client. Adjust winit window dimensions
        // on event trigger
        let closure = wasm_bindgen::closure::Closure::wrap(Box::new(move |_e: web_sys::Event| {
            let size = get_window_size();
            window.set_inner_size(size)
        }) as Box<dyn FnMut(_)>);
        client_window
            .add_event_listener_with_callback("resize", closure.as_ref().unchecked_ref())
            .unwrap();
        closure.forget();
    }

    let mut input = WinitInputHelper::new();
    let mut pixels = {
        let window_size = window.inner_size();
        let surface_texture =
            SurfaceTexture::new(window_size.width, window_size.height, window.as_ref());
        Pixels::new_async(WIDTH, HEIGHT, surface_texture)
            .await
            .expect("Pixels error")
    };

    let mut shapes = vec![
        get_cube().with_location(Pt3::new(-7.0, 2.0, 2.0)),
        get_tetrahedron().with_location(Pt3::new(-3.0, 2.0, 2.0)),
        get_icosahedron().with_location(Pt3::new(5.0, 2.0, 2.0)),
        get_suzanne().with_location(Pt3::new(7.0, 2.0, 2.0)),
        get_teapot().with_location(Pt3::new(0.0, 0.0, 2.0)),
    ];

    shapes.extend(get_floor());
    shapes.extend(get_axes());

    shapes.push(get_cube().with_scale(0.1));

    let mut world = world::World {
        world_bbox: None,
        camera: Camera {
            location: Vec3::new(0.0, 0.0, 0.5),
            rotation: Quat::from_euler_angles(0.0, 0.0, 0.0),
            rotation_angles: Vec3::zeros(),
            projection: world::Projection::Orthographic,
            fov_radians: PI / 6.0,
            target_height: 2.0,
            light_vector: -Vec3::y_axis(),
            light_vector_eulers: Vec3::zeros(),
        },
        shapes,
    };

    world.camera.recalculate_light();

    pixels.frame_mut().fill(255);

    event_loop.run(move |event, _window_target, control_flow| {
        // Draw the current frame
        if let Event::RedrawRequested(_) = event {
            let mut target = ArrayDrawTarget {
                array: pixels.frame_mut(),
                width: WIDTH as usize,
                pixels_drawn: 0,
            };
            world.rasterize(&mut target).unwrap();
            if let Err(err) = pixels.render() {
                log_error("pixels.render", err);
                *control_flow = ControlFlow::Exit;
                return;
            }
        }

        // Handle input events
        if input.update(&event) {
            // Close events
            if input.key_pressed(VirtualKeyCode::Escape) || input.close_requested() {
                *control_flow = ControlFlow::Exit;
                return;
            }

            // Camera coordinates:
            // +X -- right
            // +Y -- up
            // +Z -- out of screen

            //*control_flow = ControlFlow::Wait;
            let speed = 0.2;
            if input.key_held(VirtualKeyCode::A) {
                // left
                world.camera.location += world.camera.rotation * Vec3::new(-speed, 0.0, 0.0);
                world.camera.location.y = world.camera.target_height;
                *control_flow = ControlFlow::Poll;
            }
            if input.key_held(VirtualKeyCode::D) {
                // Right
                world.camera.location += world.camera.rotation * Vec3::new(speed, 0.0, 0.0);
                world.camera.location.y = world.camera.target_height;
                *control_flow = ControlFlow::Poll;
            }
            if input.key_held(VirtualKeyCode::W) {
                // Forward into screen
                world.camera.location += world.camera.rotation * Vec3::new(0.0, 0.0, -speed);
                world.camera.location.y = world.camera.target_height;
                *control_flow = ControlFlow::Poll;
            }
            if input.key_held(VirtualKeyCode::S) {
                // Back out of screen
                world.camera.location += world.camera.rotation * Vec3::new(0.0, 0.0, speed);
                world.camera.location.y = world.camera.target_height;
                *control_flow = ControlFlow::Poll;
            }

            if input.key_held(VirtualKeyCode::H) {
                world.camera.light_vector_eulers.x -= 0.1;
                world.camera.recalculate_light();
                *control_flow = ControlFlow::Poll;
            }
            if input.key_held(VirtualKeyCode::L) {
                world.camera.light_vector_eulers.x += 0.1;
                world.camera.recalculate_light();
                *control_flow = ControlFlow::Poll;
            }
            if input.key_held(VirtualKeyCode::J) {
                world.camera.light_vector_eulers.y -= 0.1;
                world.camera.recalculate_light();
                *control_flow = ControlFlow::Poll;
            }

            if input.key_held(VirtualKeyCode::K) {
                world.camera.light_vector_eulers.y += 0.1;
                world.camera.recalculate_light();
                *control_flow = ControlFlow::Poll;
            }

            // Compute the camera rotation.
            if input.mouse_held(0) {
                let (dx, dy) = input.mouse_diff();
                let sensitivity = 400.0;
                let dx = dx / sensitivity;
                let dx = dx as f64;
                let dy = dy / sensitivity;
                let dy = dy as f64;

                world.camera.rotation_angles += Vec3::new(dy, dx, 0.0);
                world.camera.rotation = UnitQuaternion::from_euler_angles(
                    world.camera.rotation_angles.x,
                    world.camera.rotation_angles.y,
                    world.camera.rotation_angles.z,
                );

                // let rotation = UnitQuaternion::from_euler_angles(dy, dx, 0.0);
                // world.camera.rotation *= rotation;
            }

            if input.mouse_held(2) {
                world.camera.fov_radians += f64::to_radians(input.mouse_diff().1 as f64);
            }

            if input.key_pressed(VirtualKeyCode::P) {
                world.camera.projection = match world.camera.projection {
                    world::Projection::Perspective => world::Projection::Orthographic,
                    world::Projection::Orthographic => world::Projection::Perspective,
                }
            }

            if input.key_pressed(VirtualKeyCode::C) {
                for object in world.shapes.iter_mut() {
                    object.style = match object.style {
                        world::RenderingStyle::Wireframe { line_style } => {
                            world::RenderingStyle::Unlit {
                                color: line_style.stroke_color.unwrap_or_default(),
                                line_color: None,
                            }
                        }
                        world::RenderingStyle::UnlitRandom => world::RenderingStyle::UnlitRandom,
                        world::RenderingStyle::Unlit { color, line_color } => {
                            world::RenderingStyle::ViewLight { color, line_color }
                        }
                        world::RenderingStyle::ViewLight { color, line_color } => {
                            world::RenderingStyle::GlobalLight { color, line_color }
                        }
                        world::RenderingStyle::GlobalLight {
                            color,
                            line_color: _,
                        } => world::RenderingStyle::Wireframe {
                            line_style: PrimitiveStyleBuilder::new()
                                .stroke_color(color)
                                .stroke_width(2)
                                .build(),
                        },
                    };
                }
            }

            // Resize the window
            if let Some(size) = input.window_resized() {
                if let Err(err) = pixels.resize_surface(size.width, size.height) {
                    log_error("pixels.resize_surface", err);
                    *control_flow = ControlFlow::Exit;
                    return;
                }
            }

            // Update internal state and request a redraw
            world.update();
            window.request_redraw();
        }
    });
}

fn log_error<E: std::error::Error + 'static>(method_name: &str, err: E) {
    error!("{method_name}() failed: {err}");
    for source in err.sources().skip(1) {
        error!("  Caused by: {source}");
    }
}
