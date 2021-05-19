#![allow(dead_code)]

use std::ops;
use image;

const STEP : f32 = DELTA/(N as f32);
const EPSILON : f32 = 0.0001;   // to compare floating-point values

const N : usize = 2000;   // screen is a N by N grid (middle of screen is the origin)
const DELTA : f32 = 10.;   // size of screen within the scene
const OMEGA : Pt = Pt(0., 0., 10.);   // location of viewer

const BACKGROUND : Color = Color([0., 0., 0.]);


fn main() {

    let objects : Vec<Sphere> = vec![Sphere {ctr : Pt(-3., 0., -5.), rad : 3.}, Sphere {ctr : Pt(3., 0., -5.), rad : 2.}];    // all Objects in scene
    let mut sources : Vec<Pt> = vec![Pt(0., 0., -2.), Pt(0., -5., -2.)];   // all Light Sources in scene
    let kd : Vec<Color> = vec![Color([1., 0., 0.]), Color([0., 1., 1.])];    // all Objects' kd
    let colors : Vec<Color> = vec![Color([0.5, 0.5, 0.5]), Color([0.25, 0.5, 0.25])];   // all Light Sources' color

    for i in 0..100 {
        let y = (-5.)*(i as f32 / 100.);
        sources[1] = Pt(0., y, -2.);
        let mut name = i.to_string();
        name.push_str(".png");
        run(&objects, &kd, &sources, &colors, name);
    }
}

#[derive(PartialEq, Clone, Copy, Debug)]
struct Pt(f32, f32, f32);

#[derive(PartialEq, Clone, Copy, Debug)]
struct Ve(f32, f32, f32);

#[derive(PartialEq, Clone, Copy, Debug)]
struct Ray {
    src : Pt,
    dir : Ve,
}

#[derive(PartialEq, Clone, Copy, Debug)]
struct Sphere {
    ctr : Pt,
    rad : f32,
}

#[derive(PartialEq, Clone, Copy, Debug)]
struct Color([f32;3]);

impl Ve {

    fn of(p : Pt) -> Ve {
        Ve(p.0, p.1, p.2)
    }

    fn from(a : Pt, b : Pt) -> Ve {
        Ve(b.0-a.0, b.1-a.1, b.2-a.2)
    }

    fn norm(&self) -> f32 {
        (self.dot(self)).sqrt()
    }

    fn unit(self) -> Ve {
        let self_hat = (1./self.norm()) * self;
        assert!((self_hat.norm() - 1.).abs() <= EPSILON);
        self_hat
    }

    fn dot(&self, v : &Ve) -> f32 {
        self.0*v.0 + self.1*v.1 + self.2*v.2
    }

    fn dir(a : Pt, b : Pt) -> Ve {
        (Ve::from(a, b)).unit()
    }
}

impl ops::Add<Ve> for Ve {
    type Output = Ve;

    fn add(self, u : Ve) -> Ve {
        Ve(self.0+u.0, self.1+u.1, self.2+u.2)
    }
}

impl ops::Mul<Ve> for f32 {
    type Output = Ve;

    fn mul(self, u : Ve) -> Ve {
        Ve(self * u.0, self * u.1, self * u.2)
    }
}

impl Pt {

    fn of(u : Ve) -> Pt {
        Pt(u.0, u.1, u.2)
    }

    fn visible(self, objects : &[Sphere], j : usize, src : Pt) -> bool{
    
        let mut visibility = true;
        let r = Ray::from_pts(src, self);
        let t1 = Ve::from(src, self).norm();
    
        for (i, obj) in objects.iter().enumerate() {
            if i != j {
                if let Some((_, t2)) = obj.intersect(r) {
                    if t2 < t1 {
                        visibility = false;
                        break;
                    } else {} 
                } else {}
            } else {
                if !(obj.over_horizon(self, src)) {
                    visibility = false;
                    break;
                } else {}
            }
        }
    
        visibility
    
    }

    fn diffused_color(self, objects : &[Sphere], kd : &[Color], sources : &[Pt], colors : &[Color], j : usize) -> Color {
        let mut total_color = Color([0., 0., 0.]);
        for (k, src) in sources.iter().enumerate() {
            if self.visible(objects, j, *src) {
                total_color = total_color + Ray::from_pts(*src, self).diffused_color(colors[k], Ve::dir(self, objects[j].ctr), kd[j])
            }
        }
        total_color
    }
}

impl Ray {

    fn from_pts(a : Pt, b : Pt) -> Ray {
        Ray {src : a,
             dir :Ve::dir(a, b)
            }
    }

    fn pt(&self, t : f32) -> Pt {
        assert!(t >= 0.);
        let Ray {src : s, dir : d_hat} = *self;
        Pt::of(Ve::of(s) + t*d_hat)
    }

    fn diffused_color(&self, cs : Color,  n : Ve, kd : Color) -> Color {
        // uses Lambert's law
        assert!((n.norm() - 1.).abs() <= EPSILON);
        let cos_theta = n.dot(&self.dir.unit());
        Color([cs.0[0]*kd.0[0]*(cos_theta), cs.0[1]*kd.0[1]*(cos_theta), cs.0[2]*kd.0[2]*(cos_theta)])
    }

    fn intercept(self, objects : &[Sphere]) -> Option<(Pt, usize)> {
        let mut t = f32::INFINITY;
        let mut p = Pt(0., 0., 0.);
        let mut i = 0;
        for (k, obj) in objects.iter().enumerate() {
            if let Some((p2, t2)) = obj.intersect(self) {
                if t > t2 {
                    t = t2;
                    p = p2;
                    i = k;
                } else {}
            } else {}
        }

        if t == f32::INFINITY {
            None
        } else {
            Some((p, i))
        }
    }

}

impl Sphere {
    fn from(a : Pt, b : Pt) -> Sphere{
        let r = Ve::from(a, b).norm();
        Sphere {ctr : a, rad : r}
    }

    fn intersect(self, r : Ray) -> Option<(Pt, f32)> {
        // Solves 𝑡2 +2𝑡( ⃗𝑢⋅⃗⃗⃗⃗⃗⃗⃗⃗⃗⃗⃗⃗𝐶𝐴)+‖⃗⃗⃗⃗⃗⃗⃗⃗⃗⃗⃗⃗𝐶𝐴‖2 −𝑟2 = 0 for t using the discriminant of this polynomial
        let a : f32 = 1. ;
        let b = 2. * (r.dir).dot(&Ve::from(self.ctr, r.src));
        let c = Ve::from(self.ctr, r.src).norm().powi(2) - (self.rad).powi(2);

        let discriminant = (b*b) - 4.*(a*c);

        if discriminant >= 0. {
           let  t1 = (-b - discriminant.sqrt())/2.*a;
           let t2 = (-b + discriminant.sqrt())/2.*a;
           if t2 >= t1 && t1 >= 0. {
               Some((r.pt(t1), t1))
           } else if t2 >= 0. {
                Some((r.pt(t2), t2))
           } else {
               None
           }
        } else {
            None
        }

    }

    fn over_horizon(self, p : Pt, src : Pt) -> bool {
        (Ve::from(self.ctr, p)).dot(&Ve::from(p, src)) >= 0.
    }
}

impl ops::Add<Color> for Color {
    type Output = Color;

    fn add(self, u : Color) -> Color {
        Color([self.0[0]+u.0[0], self.0[1]+u.0[1], self.0[2]+u.0[2]])
    }
}

fn grid(i : usize, j : usize) -> Pt {
    Pt(((i as f32) + 0.5 - (N/2) as f32)*STEP, ((j as f32) + 0.5 - (N/2) as f32)*STEP, 0.)
}

fn shoot_ray(i : usize, j : usize) -> Ray {
    Ray::from_pts(OMEGA, grid(i, j))
}

fn run(objects : &[Sphere], kd : &[Color], sources : &[Pt], colors : &[Color], name : String){

    assert!(kd.len() == objects.len() && colors.len() == sources.len());

    let mut imgbuf = image::ImageBuffer::new(N as u32, N as u32);

    for i in 0..N {
        for j in 0..N {
            let r = shoot_ray(i, j);
            if let Some((p, k)) = r.intercept(objects) {
                let color = p.diffused_color(objects, kd, sources, colors, k);
                let [r_flt, g_flt, b_flt] = color.0;
                assert!(r_flt>=0. && r_flt<=1.);
                assert!(g_flt>=0. && g_flt<=1.);
                assert!(b_flt>=0. && b_flt<=1.);
                let r = (r_flt*255.).round() as u8;
                let g = (g_flt*255.).round() as u8;
                let b = (b_flt*255.).round() as u8;
                imgbuf.put_pixel(i as u32, j as u32, image::Rgb([r, g, b]));
            } else {}
        }
    }
    imgbuf.save(name).unwrap();
}

//-----------------------------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    #![warn(dead_code)]
    use super::*;

    #[test]
    fn operations() {
        let u = Ve(1.5, 4., 0.);
        let v = Ve(1.5, 6., 1.75);

        assert_eq!(u+v, Ve(3., 10., 1.75));

        assert_eq!(1.5*u, Ve(2.25, 6., 0.));
    }

    #[test]
    fn ve_from_pts() {

        let a = Pt(1., 2., 3.);
        let b = Pt(-3.5, 5., 2.);

        let u = Ve::from(a, b);

        assert!(u == Ve(-4.5, 3., -1.));
    }

    #[test]
    fn dot_prod(){
        let u = Ve(2., 1.5, 0.);
        let v = Ve(1.5, 1.5, 4.);

        let x = u.dot(&v);

        assert_eq!(x, 5.25);
    }

    #[test]
    fn unit_vectors() {
        let u = Ve(3., 0., 0.);

        let u_hat = u.unit();

        assert_eq!(u_hat, Ve(1., 0., 0.));
    }

    #[test]
    fn pt_of_ray() {
        let t : f32 = 3.5;
        let r = Ray {src : Pt(1., 1., 1.), dir : Ve(1., 0., 1.)};

        let p = r.pt(t);

        assert_eq!(p, Pt(4.5, 1., 4.5));
    }

    #[test]
    fn direction() {
        let a = Pt(1., 1., 0.);
        let b = Pt(1., 1., 3.);

        let u_hat = Ve::dir(a, b);

        assert_eq!(u_hat, Ve(0., 0., 1.));
    }

    #[test]
    fn sphere_from_pts() {
        let sph = Sphere::from(Pt(0., 3., 0.), Pt(0., 4., 0.));
        assert_eq!(sph, Sphere {ctr : Pt(0., 3., 0.), rad : 1.});
    }

    #[test]
    fn ray_from_pts() {
        let r = Ray::from_pts(Pt(0., 0., 0.), Pt(3., 0., 0.));
        assert_eq!(r, Ray {src : Pt(0., 0., 0.), dir : Ve(1., 0., 0.)});
    }

    #[test]
    fn intersection_sphere() {
        let sph = Sphere::from(Pt(0., 0., 0.), Pt(1., 0., 0.));
        let r = Ray{src : Pt(2., 0., 0.), dir : Ve(-1., 0., 0.)};

        assert_eq!(sph.intersect(r), Some((Pt(1., 0., 0.), 1.)));

        let sph = Sphere::from(Pt(0., 0., 0.), Pt(1., 0., 0.));
        let r = Ray{src : Pt(0., 0., -3.), dir : Ve(0., 3., 1.).unit()};

        assert_eq!(sph.intersect(r), None);
    }

    #[test]
    fn horizon() {
        let sph = Sphere::from(Pt(0., 0., 0.), Pt(1., 0., 0.));
        let p = Pt(0., 1., 0.);

        let src1 = Pt(1., 1.5, 0.);
        let src2 = Pt(2., 0., 0.);

        assert_eq!(sph.over_horizon(p, src1), true);
        assert_eq!(sph.over_horizon(p, src2), false);
    }

    #[test]
    fn visibility() {
        let sph1 = Sphere::from(Pt(0., 0., 0.), Pt(1., 0., 0.));
        let sph2 = Sphere::from(Pt(0., 3., 0.), Pt(0., 4., 0.));

        let objects = vec![sph1, sph2];
        let j : usize = 0;
        

        let p = Pt(0., 1., 0.);

        let src1 = Pt(0., 1.5, 0.);
        let src2 = Pt(0., 5., 0.);
        
        assert_eq!(p.visible(&objects, j, src1), true);
        assert_eq!(p.visible(&objects, j, src2), false);
   }

   #[test]
   fn interceptions() {
    let sph1 = Sphere {ctr : Pt(0., 0., 0.), rad : 1.};
    let sph2 = Sphere {ctr : Pt(0., 0., 1.), rad : 1.};

    let objects = vec![sph1, sph2];

    let r1 = Ray {src : Pt(0., 0., -3.), dir : Ve(0., 0., 1.).unit()};
    let r2 = Ray {src : Pt(0., 0., -3.), dir : Ve(0., 3., 1.).unit()};

    assert_eq!(r1.intercept(&objects), Some((Pt(0., 0., -1.) , 0)));
    assert_eq!(r2.intercept(&objects), None);
   }

   #[test]
   fn screen_coords() {
       assert_eq!(grid(N/2, N/2), Pt(0.5*STEP, 0.5*STEP, 0.));
       //assert_eq!(grid(0, 0), Pt(-0.5*(DELTA as f32), 0.5*(DELTA as f32), 0.)); //seems to work : floats are close enough
   }

   #[test]
   fn pt_color() {
       let p = Pt(0., 1., 0.);
       let objects = vec![Sphere {ctr : Pt(0., 0., 0.), rad : 1.}];
       let kd = vec![Color([1., 0., 0.])];
       let sources = vec![Pt(0., 2., 0.), Pt(0., -2., 0.)];
       let colors = vec![Color([1., 1., 1.]), Color([0., 0.5, 1.])]; // Note : Do not reproduce, colors should add up to less than 1. on each channel !!!
       assert_eq!(p.diffused_color(&objects, &kd, &sources, &colors, 0), Color([1., 0., 0.]));
   }

   #[test]
   fn summing_colors() {
       let a = Color([0., 0.5, 0.5]);
       let b = Color([1., 0.5, 0.]);
       assert_eq!(a+b, Color([1., 1., 0.5]));
   }

}