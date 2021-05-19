pub fn lin_grad(c1 : [u8;3], c2 : [u8;3], n : u32) -> Vec<[u8;3]> {
    let mut grad = Vec::new();
    for i in 0..= n {
        let t = (i as f32)/(n as f32);
        let [r1, g1, b1] = c1;
        let [r2, g2, b2] = c2;
        let r = ((r1 as f32) * (1.0 - t) + (r2 as f32)*t).round() as u8;
        let g = ((g1 as f32) * (1.0 - t) + (g2 as f32)*t).round() as u8;
        let b = ((b1 as f32) * (1.0 - t) + (b2 as f32)*t).round() as u8;
        grad.push([r, g, b]);
    }
    grad
}

pub fn multi_lin_grad(colors : &[[u8;3]], n : u32) -> Vec<[u8;3]> {
    let l = colors.len();
    let n_out = ((n as f32)/(l as f32)).round() as u32;
    let mut grad = Vec::new();
    for i in 0..(l-1) {
        grad.extend(lin_grad(colors[i], colors[i+1], n_out));
    }
    while grad.len() <= n as usize {
        grad.push(colors[l-1])
    }
    grad
}
