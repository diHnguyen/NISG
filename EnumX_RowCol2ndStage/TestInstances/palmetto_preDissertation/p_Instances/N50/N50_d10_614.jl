global arcs = [1 2; 1 11; 1 21; 1 26; 1 31; 1 32; 2 12; 2 33; 2 39; 3 6; 3 17; 3 26; 3 34; 4 3; 4 16; 4 22; 4 32; 4 38; 4 41; 4 48; 5 9; 5 10; 5 31; 5 34; 5 45; 5 47; 6 12; 6 17; 6 26; 7 13; 7 15; 7 22; 7 26; 7 27; 7 33; 7 41; 7 50; 8 5; 8 18; 8 36; 8 38; 9 3; 9 30; 9 36; 9 37; 9 42; 9 44; 9 45; 9 46; 10 9; 10 14; 10 20; 10 21; 11 13; 11 43; 11 46; 12 4; 12 20; 12 22; 12 30; 12 40; 12 42; 12 44; 13 14; 13 26; 13 33; 14 42; 15 20; 15 27; 15 31; 15 34; 15 36; 15 43; 16 10; 16 12; 16 33; 16 44; 16 49; 17 3; 17 21; 17 29; 17 36; 17 46; 17 48; 18 4; 18 8; 18 23; 18 25; 18 28; 18 29; 18 38; 18 45; 19 12; 19 41; 19 45; 19 50; 20 11; 20 14; 20 23; 20 39; 21 4; 21 16; 21 20; 21 26; 21 45; 22 5; 22 47; 23 9; 23 12; 23 16; 23 24; 23 30; 24 6; 24 28; 24 46; 24 47; 25 9; 25 12; 25 16; 25 33; 25 50; 26 2; 26 13; 26 31; 26 40; 26 42; 27 9; 27 14; 27 22; 27 25; 27 28; 27 34; 27 40; 27 43; 28 29; 28 33; 28 50; 29 5; 29 15; 29 26; 29 31; 29 37; 30 8; 30 18; 30 21; 30 42; 30 44; 30 45; 31 7; 31 20; 31 30; 31 43; 32 2; 32 15; 32 33; 33 15; 33 24; 33 37; 33 47; 34 4; 34 18; 34 38; 35 10; 35 11; 35 16; 36 18; 36 24; 36 27; 36 31; 36 33; 36 42; 37 2; 38 39; 38 50; 39 3; 39 26; 39 38; 39 40; 40 2; 40 8; 40 9; 40 13; 40 16; 40 23; 40 41; 40 44; 41 2; 41 10; 42 9; 42 21; 42 29; 42 38; 42 50; 43 8; 43 16; 43 28; 43 40; 43 42; 43 47; 44 23; 44 28; 44 47; 44 49; 45 10; 45 25; 45 28; 45 35; 45 44; 46 5; 46 7; 46 29; 47 4; 47 22; 47 25; 47 29; 47 37; 48 30; 48 34; 48 47; 49 2; 49 11; 49 19; 49 26; 49 32; 49 42]
global d_x = [2.0, 6.0, 2.0, 8.0, 7.0, 9.0, 5.0, 4.0, 3.0, 4.0, 1.0, 4.0, 8.0, 9.0, 5.0, 5.0, 4.0, 5.0, 10.0, 1.0, 3.0, 6.0, 3.0, 2.0, 8.0, 8.0, 9.0, 7.0, 6.0, 3.0, 10.0, 1.0, 10.0, 2.0, 6.0, 10.0, 7.0, 9.0, 3.0, 5.0, 3.0, 4.0, 2.0, 10.0, 8.0, 2.0, 4.0, 3.0, 8.0, 4.0, 1.0, 4.0, 10.0, 10.0, 6.0, 6.0, 9.0, 5.0, 6.0, 1.0, 2.0, 7.0, 4.0, 9.0, 10.0, 2.0, 1.0, 2.0, 8.0, 5.0, 5.0, 4.0, 4.0, 8.0, 1.0, 5.0, 6.0, 1.0, 1.0, 8.0, 6.0, 5.0, 1.0, 1.0, 2.0, 3.0, 9.0, 8.0, 1.0, 5.0, 6.0, 2.0, 2.0, 4.0, 8.0, 9.0, 6.0, 3.0, 1.0, 8.0, 7.0, 7.0, 2.0, 4.0, 8.0, 4.0, 7.0, 1.0, 6.0, 5.0, 7.0, 1.0, 6.0, 7.0, 1.0, 9.0, 5.0, 7.0, 6.0, 7.0, 10.0, 8.0, 1.0, 4.0, 1.0, 9.0, 10.0, 9.0, 9.0, 7.0, 5.0, 10.0, 6.0, 6.0, 7.0, 9.0, 9.0, 7.0, 1.0, 8.0, 7.0, 9.0, 10.0, 2.0, 8.0, 4.0, 3.0, 1.0, 1.0, 9.0, 8.0, 4.0, 7.0, 9.0, 3.0, 7.0, 1.0, 6.0, 2.0, 4.0, 3.0, 1.0, 2.0, 1.0, 10.0, 9.0, 8.0, 2.0, 6.0, 9.0, 3.0, 8.0, 3.0, 4.0, 9.0, 2.0, 8.0, 3.0, 1.0, 2.0, 5.0, 8.0, 3.0, 10.0, 9.0, 7.0, 6.0, 6.0, 3.0, 5.0, 10.0, 6.0, 3.0, 3.0, 4.0, 10.0, 3.0, 4.0, 4.0, 3.0, 3.0, 1.0, 2.0, 5.0, 9.0, 5.0, 10.0, 3.0, 2.0, 1.0, 6.0, 4.0, 3.0, 4.0, 6.0, 10.0, 2.0, 4.0, 7.0, 10.0, 10.0, 1.0, 7.0, 4.0, 6.0]
global b_x = 5
global d_y = [8.0, 7.0, 10.0, 7.0, 10.0, 9.0, 2.0, 6.0, 1.0, 10.0, 5.0, 4.0, 4.0, 3.0, 5.0, 9.0, 9.0, 2.0, 6.0, 5.0, 9.0, 4.0, 8.0, 1.0, 3.0, 10.0, 4.0, 8.0, 10.0, 7.0, 5.0, 4.0, 5.0, 4.0, 6.0, 8.0, 1.0, 10.0, 8.0, 1.0, 6.0, 7.0, 6.0, 10.0, 6.0, 9.0, 9.0, 9.0, 3.0, 4.0, 8.0, 9.0, 4.0, 6.0, 8.0, 7.0, 5.0, 8.0, 10.0, 4.0, 1.0, 9.0, 6.0, 8.0, 4.0, 2.0, 9.0, 9.0, 9.0, 10.0, 6.0, 3.0, 6.0, 2.0, 2.0, 10.0, 1.0, 9.0, 3.0, 8.0, 1.0, 6.0, 5.0, 1.0, 7.0, 10.0, 2.0, 8.0, 9.0, 5.0, 9.0, 4.0, 4.0, 4.0, 7.0, 2.0, 10.0, 2.0, 7.0, 3.0, 4.0, 7.0, 4.0, 5.0, 1.0, 8.0, 6.0, 4.0, 8.0, 4.0, 8.0, 3.0, 8.0, 8.0, 4.0, 1.0, 5.0, 8.0, 4.0, 9.0, 9.0, 5.0, 2.0, 5.0, 8.0, 9.0, 7.0, 4.0, 2.0, 7.0, 9.0, 4.0, 10.0, 4.0, 9.0, 6.0, 3.0, 10.0, 9.0, 8.0, 10.0, 10.0, 3.0, 4.0, 5.0, 1.0, 3.0, 9.0, 1.0, 7.0, 1.0, 10.0, 9.0, 1.0, 2.0, 5.0, 7.0, 6.0, 1.0, 2.0, 9.0, 5.0, 7.0, 3.0, 8.0, 2.0, 10.0, 9.0, 10.0, 5.0, 3.0, 4.0, 3.0, 5.0, 7.0, 4.0, 4.0, 6.0, 9.0, 6.0, 7.0, 5.0, 4.0, 7.0, 7.0, 8.0, 7.0, 7.0, 9.0, 2.0, 2.0, 5.0, 7.0, 3.0, 6.0, 10.0, 5.0, 10.0, 10.0, 9.0, 1.0, 10.0, 1.0, 10.0, 5.0, 3.0, 1.0, 6.0, 3.0, 5.0, 8.0, 8.0, 6.0, 9.0, 8.0, 7.0, 10.0, 3.0, 5.0, 8.0, 4.0, 1.0, 2.0, 3.0, 10.0]
global b_y = 10
global p = [0.043, 0.852, 0.855, 0.417, 0.872, 0.602, 0.803, 0.905, 0.209, 0.536, 0.892, 0.941, 0.688, 0.679, 0.442, 0.615, 0.882, 0.066, 0.334, 0.623, 0.205, 0.771, 0.511, 0.112, 0.864, 0.722, 0.649, 0.089, 0.107, 0.421, 0.557, 0.626, 0.825, 0.634, 0.425, 0.43, 0.918, 0.619, 0.32, 0.638, 0.453, 0.649, 0.37, 0.78, 0.809, 0.549, 0.202, 0.642, 0.053, 0.315, 0.314, 0.958, 0.084, 0.237, 0.84, 0.335, 0.243, 0.875, 0.612, 0.993, 0.499, 0.481, 0.085, 0.799, 0.54, 0.584, 0.769, 0.106, 0.009, 0.231, 0.374, 0.133, 0.872, 0.155, 0.911, 0.27, 0.599, 0.465, 0.692, 0.828, 0.596, 0.818, 0.48, 0.635, 0.86, 0.442, 0.606, 0.065, 0.44, 0.149, 0.556, 0.028, 0.317, 0.137, 0.983, 0.859, 0.733, 0.766, 0.126, 0.567, 0.229, 0.775, 0.982, 0.846, 0.199, 0.298, 0.75, 0.483, 0.029, 0.104, 0.995, 0.624, 0.854, 0.927, 0.617, 0.167, 0.913, 0.102, 0.959, 0.582, 0.215, 0.633, 0.685, 0.599, 0.531, 0.179, 0.55, 0.447, 0.19, 0.374, 0.014, 0.978, 0.636, 0.519, 0.261, 0.658, 0.986, 0.776, 0.541, 0.89, 0.439, 0.739, 0.357, 0.048, 0.203, 0.269, 0.598, 0.842, 0.168, 0.724, 0.49, 0.448, 0.066, 0.956, 0.113, 0.647, 0.746, 0.091, 0.804, 0.163, 0.988, 0.426, 0.867, 0.72, 0.753, 0.76, 0.025, 0.508, 0.498, 0.278, 0.062, 0.543, 0.105, 0.4, 0.867, 0.022, 0.577, 0.194, 0.49, 0.543, 0.315, 0.768, 0.971, 0.051, 0.762, 0.892, 0.713, 0.474, 0.485, 0.668, 0.937, 0.777, 0.861, 0.432, 0.085, 0.028, 0.274, 0.642, 0.397, 0.239, 0.811, 0.519, 0.2, 0.53, 0.317, 0.795, 0.854, 0.029, 0.239, 0.281, 0.674, 0.615, 0.986, 0.647, 0.665, 0.948, 0.632, 0.858, 0.383, 0.61, 0.859, 0.802, 0.536, 0.7, 0.393]
global q = [0.202, 0.996, 0.984, 0.895, 0.942, 0.863, 0.993, 0.905, 0.898, 0.685, 0.982, 0.993, 0.928, 0.82, 0.904, 0.987, 0.95, 0.783, 0.975, 0.636, 0.726, 0.857, 0.565, 0.987, 0.908, 0.864, 0.784, 0.591, 0.943, 0.848, 0.76, 0.997, 0.836, 0.655, 0.442, 0.976, 0.992, 0.982, 0.55, 0.695, 0.504, 0.881, 0.433, 0.862, 0.85, 0.662, 0.66, 0.817, 0.097, 0.906, 0.709, 0.998, 0.462, 0.809, 0.858, 0.73, 0.834, 0.908, 0.887, 0.998, 0.874, 0.939, 0.135, 0.842, 0.587, 0.874, 0.834, 0.322, 0.895, 0.534, 0.523, 0.502, 0.93, 0.22, 0.986, 0.551, 0.697, 0.868, 0.723, 0.921, 0.613, 0.945, 0.714, 0.887, 0.986, 0.813, 0.793, 0.148, 0.889, 0.422, 0.652, 0.31, 0.804, 0.744, 0.991, 0.899, 0.924, 0.897, 0.319, 0.889, 0.984, 0.862, 0.992, 0.85, 0.841, 0.583, 0.87, 0.556, 0.03, 0.536, 0.999, 0.926, 0.958, 0.964, 0.918, 0.684, 0.916, 0.448, 0.985, 0.924, 0.522, 0.796, 0.767, 0.707, 0.569, 0.848, 0.951, 0.936, 0.422, 0.947, 0.914, 0.986, 0.712, 0.922, 0.792, 0.718, 0.991, 0.915, 0.845, 0.993, 0.611, 0.859, 0.389, 0.959, 0.594, 0.706, 0.937, 0.913, 0.907, 0.762, 0.932, 0.924, 0.068, 0.977, 0.796, 0.789, 0.896, 0.97, 0.914, 0.234, 0.99, 0.695, 0.881, 0.761, 0.913, 0.995, 0.049, 0.659, 0.575, 0.5, 0.315, 0.697, 0.844, 0.444, 0.91, 0.54, 0.657, 0.323, 0.727, 0.74, 0.847, 0.803, 0.977, 0.94, 0.852, 0.99, 0.766, 0.664, 0.789, 0.672, 0.942, 0.827, 0.937, 0.853, 0.912, 0.633, 0.77, 0.883, 0.85, 0.413, 0.977, 0.781, 0.724, 0.555, 0.328, 0.817, 0.927, 0.431, 0.493, 0.922, 0.798, 0.618, 0.992, 0.648, 0.833, 0.989, 0.745, 0.871, 0.639, 0.844, 0.889, 0.928, 0.77, 0.808, 0.937]
global origin = 1
global destination = 50