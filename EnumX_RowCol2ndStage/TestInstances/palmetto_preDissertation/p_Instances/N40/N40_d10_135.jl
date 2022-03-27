global arcs = [1 10; 1 19; 1 36; 2 9; 2 11; 2 20; 2 35; 2 38; 3 4; 3 6; 3 8; 3 13; 3 34; 4 30; 5 9; 5 19; 5 38; 6 18; 6 21; 6 40; 7 5; 7 13; 8 6; 8 7; 8 23; 8 33; 9 5; 9 6; 9 23; 9 34; 10 26; 10 34; 11 7; 11 9; 11 28; 11 39; 12 2; 12 26; 12 31; 12 35; 12 40; 13 4; 13 22; 13 29; 13 31; 14 26; 14 39; 15 5; 15 12; 15 21; 15 24; 15 29; 16 11; 16 17; 16 32; 16 34; 16 36; 17 3; 17 7; 17 19; 17 22; 17 29; 17 30; 17 31; 18 13; 18 29; 18 36; 18 40; 19 17; 19 35; 20 27; 20 33; 21 11; 21 20; 21 24; 21 25; 21 37; 21 40; 22 8; 22 11; 22 16; 22 23; 22 24; 22 38; 23 4; 23 25; 23 28; 23 34; 24 7; 24 32; 25 3; 25 5; 25 13; 25 15; 25 40; 26 16; 26 39; 27 8; 27 26; 27 30; 27 31; 27 40; 28 9; 28 21; 28 30; 28 34; 28 36; 28 38; 28 39; 29 19; 29 35; 29 37; 30 33; 31 6; 31 7; 31 13; 31 14; 31 20; 31 35; 31 37; 32 4; 32 18; 32 25; 33 4; 33 10; 33 18; 33 25; 34 9; 34 11; 34 13; 35 36; 36 11; 36 20; 36 28; 36 30; 36 39; 37 5; 37 8; 37 18; 38 12; 38 32; 39 2; 39 6; 39 8; 39 10; 39 13; 39 32; 39 33; 39 34]
global d_x = [1.0, 6.0, 2.0, 7.0, 3.0, 4.0, 9.0, 4.0, 10.0, 2.0, 4.0, 7.0, 2.0, 9.0, 5.0, 3.0, 8.0, 2.0, 5.0, 8.0, 5.0, 1.0, 7.0, 10.0, 2.0, 1.0, 3.0, 6.0, 6.0, 3.0, 7.0, 10.0, 8.0, 6.0, 6.0, 8.0, 10.0, 9.0, 3.0, 7.0, 10.0, 6.0, 1.0, 10.0, 2.0, 6.0, 7.0, 7.0, 9.0, 2.0, 10.0, 9.0, 7.0, 5.0, 4.0, 8.0, 4.0, 2.0, 3.0, 8.0, 9.0, 1.0, 8.0, 3.0, 8.0, 9.0, 10.0, 6.0, 9.0, 1.0, 10.0, 8.0, 9.0, 6.0, 10.0, 5.0, 10.0, 5.0, 6.0, 8.0, 8.0, 2.0, 9.0, 3.0, 4.0, 9.0, 8.0, 6.0, 4.0, 3.0, 5.0, 8.0, 2.0, 5.0, 8.0, 8.0, 7.0, 5.0, 1.0, 6.0, 9.0, 9.0, 5.0, 7.0, 5.0, 5.0, 7.0, 1.0, 5.0, 9.0, 2.0, 7.0, 2.0, 8.0, 10.0, 2.0, 9.0, 5.0, 2.0, 8.0, 1.0, 1.0, 5.0, 3.0, 2.0, 7.0, 1.0, 4.0, 5.0, 1.0, 4.0, 8.0, 5.0, 8.0, 3.0, 5.0, 5.0, 8.0, 3.0, 10.0, 3.0, 9.0, 2.0, 10.0, 6.0, 3.0, 2.0, 9.0, 7.0]
global b_x = 5
global d_y = [7.0, 3.0, 10.0, 8.0, 2.0, 4.0, 4.0, 2.0, 4.0, 3.0, 8.0, 9.0, 4.0, 6.0, 5.0, 3.0, 8.0, 3.0, 2.0, 7.0, 3.0, 7.0, 4.0, 3.0, 9.0, 6.0, 5.0, 3.0, 8.0, 10.0, 10.0, 8.0, 8.0, 7.0, 9.0, 2.0, 6.0, 10.0, 4.0, 8.0, 4.0, 7.0, 5.0, 10.0, 7.0, 4.0, 2.0, 3.0, 2.0, 1.0, 9.0, 9.0, 10.0, 3.0, 9.0, 3.0, 5.0, 5.0, 6.0, 10.0, 3.0, 10.0, 8.0, 4.0, 8.0, 9.0, 8.0, 1.0, 9.0, 4.0, 2.0, 9.0, 8.0, 1.0, 8.0, 10.0, 6.0, 9.0, 7.0, 2.0, 8.0, 4.0, 10.0, 3.0, 1.0, 4.0, 5.0, 3.0, 7.0, 3.0, 9.0, 6.0, 1.0, 3.0, 1.0, 5.0, 6.0, 6.0, 5.0, 6.0, 9.0, 3.0, 1.0, 5.0, 7.0, 7.0, 5.0, 10.0, 7.0, 8.0, 6.0, 4.0, 6.0, 10.0, 10.0, 4.0, 2.0, 3.0, 5.0, 2.0, 2.0, 7.0, 3.0, 7.0, 3.0, 2.0, 5.0, 4.0, 6.0, 8.0, 4.0, 5.0, 8.0, 2.0, 10.0, 10.0, 7.0, 2.0, 1.0, 9.0, 1.0, 10.0, 1.0, 4.0, 10.0, 3.0, 8.0, 9.0, 5.0]
global b_y = 10
global p = [0.658, 0.032, 0.15, 0.632, 0.493, 0.536, 0.585, 0.669, 0.779, 0.837, 0.359, 0.305, 0.577, 0.474, 0.032, 0.248, 0.904, 0.604, 0.161, 0.52, 0.819, 0.096, 0.643, 0.315, 0.03, 0.153, 0.412, 0.256, 0.99, 0.307, 0.503, 0.295, 0.619, 0.479, 0.514, 0.56, 0.324, 0.702, 0.448, 0.624, 0.845, 0.533, 0.212, 0.054, 0.691, 0.081, 0.062, 0.826, 0.233, 0.594, 0.884, 0.834, 0.623, 0.412, 0.66, 0.511, 0.843, 0.571, 0.524, 0.806, 0.957, 0.123, 0.927, 0.703, 0.184, 0.296, 0.022, 0.267, 0.884, 0.759, 0.726, 0.927, 0.837, 0.665, 0.816, 0.638, 0.372, 0.036, 0.378, 0.348, 0.91, 0.591, 0.333, 0.325, 0.16, 0.264, 0.034, 0.257, 0.461, 0.544, 0.638, 0.759, 0.798, 0.765, 0.645, 0.03, 0.217, 0.788, 0.812, 0.212, 0.439, 0.008, 0.457, 0.232, 0.327, 0.982, 0.256, 0.594, 0.006, 0.297, 0.52, 0.178, 0.411, 0.036, 0.9, 0.587, 0.221, 0.235, 0.856, 0.743, 0.485, 0.929, 0.519, 0.14, 0.395, 0.926, 0.022, 0.269, 0.341, 0.933, 0.849, 0.133, 0.323, 0.185, 0.3, 0.214, 0.454, 0.162, 0.402, 0.216, 0.028, 0.533, 0.253, 0.704, 0.512, 0.872, 0.737, 0.524, 0.725]
global q = [0.724, 0.831, 0.449, 0.721, 0.717, 0.986, 0.837, 0.97, 0.944, 0.983, 0.505, 0.586, 0.655, 0.788, 0.401, 0.299, 0.946, 0.995, 0.656, 0.779, 0.959, 0.292, 0.772, 0.767, 0.111, 0.947, 0.708, 0.325, 0.993, 0.863, 0.874, 0.858, 0.741, 0.72, 0.743, 0.663, 0.494, 0.951, 0.694, 0.888, 0.993, 0.798, 0.973, 0.147, 0.771, 0.505, 0.302, 0.885, 0.291, 0.738, 0.925, 0.976, 0.917, 0.513, 0.741, 0.587, 0.937, 0.646, 0.932, 0.987, 0.971, 0.456, 0.968, 0.902, 0.734, 0.71, 0.442, 0.72, 0.967, 0.981, 0.859, 0.981, 0.993, 0.751, 0.871, 0.839, 0.661, 0.981, 0.867, 0.456, 0.97, 0.773, 0.486, 0.445, 0.523, 0.582, 0.288, 0.729, 0.633, 0.746, 0.854, 0.786, 0.915, 0.838, 0.97, 0.84, 0.529, 0.964, 0.941, 0.827, 0.638, 0.572, 0.907, 0.845, 0.353, 0.995, 0.652, 0.738, 0.258, 0.887, 0.907, 0.503, 0.991, 0.298, 0.938, 0.858, 0.936, 0.974, 0.959, 0.892, 0.501, 0.99, 0.642, 0.283, 0.642, 0.935, 0.615, 0.27, 0.785, 0.945, 0.873, 0.848, 0.843, 0.589, 0.423, 0.316, 0.622, 0.987, 0.938, 0.528, 0.874, 0.707, 0.77, 0.915, 0.841, 0.877, 0.916, 0.934, 0.781]
global origin = 1
global destination = 40