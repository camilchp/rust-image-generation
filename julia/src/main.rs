use image;
use num_complex;
use colorscale;
use std::fs::File;

use std::env;

fn main() {

    let mut param: Vec<String> = env::args().collect();

    if param[1] == "animate" {

        let file_out = File::create("animation.gif").unwrap();
        let mut encoder = image::gif::GifEncoder::new(file_out);

        let r : f32 = param[2].parse().unwrap();
        let i : f32 = param[3].parse().unwrap();
        let c : num_complex::Complex<f32> = num_complex::Complex::new(r, i);
        let size : u32 = param[4].parse().unwrap();
        let iterations : u32 = 300;
        let mut colormap = Vec::new();
            colormap.push([0, 0, 80]);
            colormap.push([180, 255, 100]);
            colormap.push([255, 50, 0]);
            colormap.push([0, 0, 0]);
        let save = Option::<String>::None;
        let nb_frames : u32 = param[5].parse().unwrap();

        for k in 0..nb_frames {
            let a = (2.0*std::f32::consts::PI/(nb_frames as f32)) * k as f32;
            let c : num_complex::Complex<f32> = num_complex::Complex::new(0.0, a).exp() * c;
    
            let frame = image::Frame::new(julia(c, size, iterations, &colormap, &save).unwrap());
            encoder.encode_frame(frame).unwrap();
        }
    } else if param[1] == "frames" {

        let r : f32 = param[2].parse().unwrap();
        let i : f32 = param[3].parse().unwrap();
        let c : num_complex::Complex<f32> = num_complex::Complex::new(r, i);
        let size : u32 = param[4].parse().unwrap();
        let iterations : u32 = 300;
        let mut colormap = Vec::new();
            colormap.push([0, 0, 80]);
            colormap.push([180, 255, 100]);
            colormap.push([255, 50, 0]);
            colormap.push([0, 0, 0]);
        let nb_frames : u32 = param[5].parse().unwrap();
        let path = param.remove(6);
        let format = param.remove(6);

        for k in 0..nb_frames {
            let mut name = String::from(&path);
                name.push_str(&(k.to_string()));
                name.push('.');
                name.push_str(&format);
            let save = Some(name);
            let a = (2.0*std::f32::consts::PI/(nb_frames as f32)) * k as f32;
            let c : num_complex::Complex<f32> = num_complex::Complex::new(0.0, a).exp() * c;
    
            julia(c, size, iterations, &colormap, &save);
        }

    } else {
        let r : f32 = param[1].parse().unwrap();
        let i : f32 = param[2].parse().unwrap();
        let c : num_complex::Complex<f32> = num_complex::Complex::new(r, i);
        let size : u32 = param[3].parse().unwrap();
        let iterations : u32 = param[4].parse().unwrap();
        let mut colormap = Vec::new();
            colormap.push([0, 0, 0]);
            colormap.push([0, 255, 0]);
            colormap.push([0, 255, 0]);
            colormap.push([255, 255, 255]);
        let save = Some(param.remove(5));

        julia(c, size, iterations, &colormap, &save);
    }
}

fn julia(c : num_complex::Complex<f32>, size : u32, iterations : u32,  colormap : &[[u8;3]], save : &Option<String>) -> Option<image::RgbaImage> {

    let n : u32 = iterations;

    let imgx : u32 = size;
    let imgy : u32  = imgx;

    let centerx : f32 = 1.5;
    let centery : f32 = 1.5;

    let scalex : f32 = 3.0/ imgx as f32;
    let scaley : f32 = 3.0/ imgy as f32;

    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);

    let grad = colorscale::multi_lin_grad(colormap, n);

    for x in 0..imgx {
        for y in 0..imgy {

            let cx = x as f32 *scalex - centerx;
            let cy = y as f32 *scaley - centery;

            let mut z = num_complex::Complex::new(cx, cy);

            let mut i = 0;
            while i < n && z.norm() <= 2.0 {
                z = z*z + c;
                i += 1;
            }

            let pixel = imgbuf.get_pixel_mut(x, y);

            *pixel = image::Rgb(grad[(i as usize)]);
        }
    }

    if let Some(name) = save {
        imgbuf.save(name).unwrap();
        None
    } else {
        Some(image::DynamicImage::ImageRgb8(imgbuf).to_rgba8())
    }
}