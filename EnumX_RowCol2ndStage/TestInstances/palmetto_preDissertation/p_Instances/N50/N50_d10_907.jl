global arcs = [1 14; 1 15; 1 19; 1 26; 1 39; 2 13; 2 24; 2 28; 2 39; 2 47; 3 24; 3 33; 3 35; 3 40; 3 47; 4 14; 4 26; 4 31; 4 35; 4 36; 4 45; 5 3; 5 18; 5 20; 5 23; 5 27; 5 35; 6 20; 6 26; 6 30; 6 48; 7 2; 7 11; 7 18; 7 24; 7 41; 7 43; 8 22; 8 26; 8 35; 8 38; 8 42; 8 43; 9 11; 9 17; 9 24; 9 33; 9 47; 10 4; 10 17; 11 38; 12 13; 12 27; 12 28; 12 35; 13 6; 13 12; 13 19; 13 36; 13 40; 13 41; 13 42; 13 44; 13 50; 14 13; 14 28; 15 2; 15 4; 15 23; 15 31; 15 34; 15 41; 16 5; 16 7; 16 11; 16 12; 16 27; 16 40; 16 41; 17 4; 17 6; 17 12; 17 15; 17 21; 17 26; 17 49; 18 2; 18 10; 18 19; 18 24; 19 4; 19 5; 19 9; 19 16; 19 30; 19 46; 19 48; 20 2; 20 24; 20 47; 21 3; 21 7; 21 24; 21 26; 21 32; 21 39; 21 42; 21 44; 21 50; 22 5; 22 10; 22 15; 22 23; 22 32; 22 36; 22 42; 22 47; 23 9; 23 21; 23 24; 23 41; 23 42; 23 45; 23 47; 24 9; 24 42; 24 47; 25 16; 25 17; 25 27; 25 43; 25 48; 26 10; 26 15; 26 18; 26 20; 26 28; 26 29; 26 30; 26 40; 26 48; 27 8; 27 29; 27 32; 27 41; 27 44; 28 4; 28 7; 28 8; 28 12; 28 13; 28 15; 28 18; 28 19; 28 26; 29 2; 29 23; 29 36; 29 39; 30 3; 30 11; 30 15; 30 20; 30 22; 30 40; 30 48; 31 10; 31 17; 31 27; 31 30; 31 35; 32 3; 32 20; 32 45; 32 46; 32 50; 33 4; 33 34; 33 39; 33 44; 34 6; 34 7; 34 23; 34 44; 34 50; 35 24; 35 32; 35 37; 35 39; 35 40; 35 44; 35 45; 36 7; 36 17; 36 21; 36 22; 36 23; 36 24; 36 25; 36 42; 37 26; 37 28; 37 30; 37 35; 38 2; 38 3; 38 9; 38 12; 38 18; 38 33; 38 34; 38 36; 38 47; 39 6; 39 10; 39 17; 39 19; 39 20; 39 28; 39 30; 39 40; 40 10; 40 20; 41 14; 41 15; 41 22; 41 25; 41 35; 41 37; 41 49; 42 7; 42 26; 42 41; 42 48; 43 4; 43 11; 43 17; 43 20; 43 45; 43 48; 44 9; 44 10; 44 19; 44 21; 44 31; 44 38; 44 46; 45 39; 45 43; 45 49; 46 2; 46 5; 46 19; 46 23; 46 29; 46 35; 47 18; 47 20; 48 4; 48 36; 48 45; 48 49; 49 7; 49 15; 49 38; 49 46]
global d_x = [3.0, 5.0, 9.0, 4.0, 10.0, 5.0, 8.0, 2.0, 1.0, 9.0, 1.0, 4.0, 3.0, 10.0, 2.0, 2.0, 1.0, 7.0, 8.0, 8.0, 5.0, 1.0, 10.0, 2.0, 4.0, 2.0, 6.0, 9.0, 2.0, 4.0, 4.0, 10.0, 10.0, 1.0, 7.0, 9.0, 5.0, 2.0, 2.0, 4.0, 1.0, 7.0, 3.0, 4.0, 7.0, 5.0, 7.0, 1.0, 6.0, 5.0, 9.0, 1.0, 5.0, 8.0, 4.0, 3.0, 1.0, 7.0, 5.0, 7.0, 4.0, 8.0, 5.0, 8.0, 1.0, 2.0, 9.0, 4.0, 9.0, 9.0, 3.0, 3.0, 8.0, 8.0, 6.0, 6.0, 8.0, 9.0, 6.0, 7.0, 8.0, 1.0, 2.0, 9.0, 1.0, 5.0, 7.0, 6.0, 5.0, 4.0, 2.0, 7.0, 2.0, 8.0, 3.0, 6.0, 4.0, 7.0, 7.0, 10.0, 1.0, 9.0, 5.0, 2.0, 9.0, 2.0, 1.0, 8.0, 2.0, 9.0, 8.0, 2.0, 7.0, 4.0, 3.0, 10.0, 8.0, 6.0, 5.0, 7.0, 2.0, 10.0, 10.0, 6.0, 5.0, 10.0, 6.0, 1.0, 1.0, 6.0, 2.0, 4.0, 3.0, 5.0, 2.0, 5.0, 9.0, 7.0, 6.0, 2.0, 2.0, 8.0, 8.0, 3.0, 1.0, 8.0, 1.0, 6.0, 7.0, 2.0, 7.0, 9.0, 2.0, 1.0, 6.0, 10.0, 7.0, 8.0, 9.0, 4.0, 5.0, 1.0, 7.0, 8.0, 10.0, 8.0, 10.0, 9.0, 4.0, 2.0, 7.0, 5.0, 4.0, 8.0, 5.0, 8.0, 9.0, 8.0, 2.0, 2.0, 4.0, 5.0, 9.0, 6.0, 2.0, 6.0, 10.0, 1.0, 1.0, 10.0, 5.0, 4.0, 3.0, 10.0, 7.0, 2.0, 7.0, 5.0, 6.0, 7.0, 6.0, 8.0, 2.0, 1.0, 2.0, 10.0, 1.0, 6.0, 3.0, 10.0, 4.0, 1.0, 4.0, 2.0, 6.0, 8.0, 2.0, 5.0, 1.0, 7.0, 10.0, 1.0, 10.0, 6.0, 2.0, 10.0, 9.0, 7.0, 1.0, 8.0, 5.0, 10.0, 5.0, 1.0, 10.0, 2.0, 6.0, 5.0, 6.0, 4.0, 9.0, 10.0, 6.0, 10.0, 2.0, 4.0, 10.0, 7.0, 5.0, 2.0, 4.0, 4.0, 1.0, 2.0, 5.0, 2.0, 7.0, 1.0, 1.0, 4.0, 5.0, 2.0, 10.0, 5.0, 1.0, 1.0]
global b_x = 5
global d_y = [8.0, 9.0, 5.0, 8.0, 1.0, 6.0, 9.0, 10.0, 2.0, 1.0, 5.0, 8.0, 7.0, 1.0, 4.0, 1.0, 7.0, 1.0, 6.0, 8.0, 9.0, 5.0, 10.0, 4.0, 1.0, 10.0, 9.0, 4.0, 8.0, 6.0, 9.0, 4.0, 4.0, 8.0, 8.0, 4.0, 9.0, 6.0, 8.0, 6.0, 1.0, 8.0, 9.0, 3.0, 10.0, 9.0, 2.0, 10.0, 4.0, 9.0, 9.0, 8.0, 5.0, 1.0, 10.0, 6.0, 7.0, 3.0, 9.0, 1.0, 7.0, 4.0, 6.0, 3.0, 7.0, 2.0, 7.0, 5.0, 3.0, 1.0, 2.0, 9.0, 9.0, 10.0, 4.0, 1.0, 10.0, 3.0, 1.0, 6.0, 6.0, 3.0, 5.0, 10.0, 1.0, 4.0, 3.0, 6.0, 7.0, 5.0, 3.0, 8.0, 2.0, 1.0, 3.0, 3.0, 7.0, 2.0, 7.0, 5.0, 4.0, 1.0, 3.0, 2.0, 1.0, 10.0, 3.0, 3.0, 3.0, 2.0, 5.0, 10.0, 1.0, 10.0, 5.0, 3.0, 6.0, 2.0, 9.0, 5.0, 4.0, 9.0, 3.0, 8.0, 4.0, 9.0, 3.0, 2.0, 10.0, 9.0, 8.0, 6.0, 4.0, 9.0, 5.0, 8.0, 8.0, 8.0, 3.0, 2.0, 6.0, 1.0, 5.0, 2.0, 3.0, 6.0, 4.0, 3.0, 4.0, 1.0, 9.0, 4.0, 7.0, 9.0, 7.0, 2.0, 2.0, 5.0, 3.0, 8.0, 1.0, 10.0, 5.0, 2.0, 6.0, 4.0, 2.0, 10.0, 1.0, 5.0, 7.0, 1.0, 3.0, 7.0, 8.0, 6.0, 9.0, 10.0, 4.0, 2.0, 7.0, 3.0, 5.0, 6.0, 3.0, 10.0, 1.0, 8.0, 4.0, 6.0, 5.0, 2.0, 4.0, 5.0, 4.0, 2.0, 8.0, 7.0, 3.0, 6.0, 5.0, 5.0, 5.0, 4.0, 2.0, 2.0, 8.0, 8.0, 4.0, 1.0, 1.0, 1.0, 5.0, 5.0, 8.0, 6.0, 6.0, 4.0, 9.0, 4.0, 7.0, 3.0, 8.0, 8.0, 1.0, 9.0, 5.0, 10.0, 4.0, 2.0, 4.0, 9.0, 3.0, 9.0, 7.0, 6.0, 5.0, 5.0, 4.0, 7.0, 9.0, 4.0, 7.0, 1.0, 6.0, 5.0, 5.0, 9.0, 1.0, 10.0, 7.0, 3.0, 6.0, 4.0, 10.0, 8.0, 5.0, 7.0, 3.0, 7.0, 10.0, 10.0, 1.0, 7.0, 10.0, 6.0]
global b_y = 10
global p = [0.115, 0.503, 0.855, 0.089, 0.247, 0.904, 0.536, 0.212, 0.521, 0.037, 0.515, 0.866, 0.19, 0.941, 0.037, 0.226, 0.445, 0.276, 0.379, 0.4, 0.827, 0.775, 0.348, 0.45, 0.199, 0.719, 0.799, 0.031, 0.673, 0.646, 0.061, 0.816, 0.333, 0.606, 0.591, 0.325, 0.139, 0.164, 0.278, 0.985, 0.305, 0.807, 0.423, 0.445, 0.857, 0.37, 0.132, 0.722, 0.056, 0.542, 0.11, 0.158, 0.589, 0.67, 0.586, 0.462, 0.46, 0.071, 0.322, 0.726, 0.585, 0.113, 0.813, 0.101, 0.295, 0.958, 0.914, 0.573, 0.923, 0.681, 0.58, 0.615, 0.76, 0.75, 0.127, 0.077, 0.748, 0.137, 0.946, 0.736, 0.018, 0.563, 0.684, 0.415, 0.076, 0.052, 0.174, 0.873, 0.79, 0.389, 0.461, 0.354, 0.472, 0.72, 0.326, 0.974, 0.129, 0.736, 0.786, 0.341, 0.751, 0.91, 0.635, 0.882, 0.063, 0.075, 0.262, 0.748, 0.777, 0.145, 0.239, 0.888, 0.494, 0.628, 0.418, 0.63, 0.717, 0.971, 0.382, 0.304, 0.357, 0.59, 0.915, 0.998, 0.852, 0.755, 0.936, 0.751, 0.323, 0.515, 0.618, 0.288, 0.787, 0.656, 0.927, 0.722, 0.43, 0.953, 0.607, 0.116, 0.895, 0.378, 0.988, 0.035, 0.823, 0.353, 0.113, 0.166, 0.594, 0.81, 0.591, 0.85, 0.338, 0.77, 0.723, 0.038, 0.622, 0.867, 0.789, 0.714, 0.934, 0.632, 0.477, 0.507, 0.168, 0.704, 0.56, 0.838, 0.438, 0.876, 0.577, 0.244, 0.612, 0.542, 0.913, 0.211, 0.209, 0.571, 0.903, 0.804, 0.954, 0.462, 0.245, 0.803, 0.156, 0.015, 0.131, 0.059, 0.304, 0.018, 0.548, 0.851, 0.362, 0.582, 0.723, 0.982, 0.658, 0.985, 0.558, 0.677, 0.409, 0.643, 0.659, 0.745, 0.668, 0.005, 0.187, 0.426, 0.974, 0.399, 0.167, 0.916, 0.434, 0.82, 0.42, 0.914, 0.666, 0.155, 0.286, 0.67, 0.233, 0.89, 0.926, 0.673, 0.749, 0.649, 0.875, 0.416, 0.007, 0.567, 0.575, 0.694, 0.387, 0.958, 0.583, 0.368, 0.295, 0.248, 0.75, 0.65, 0.22, 0.373, 0.379, 0.23, 0.657, 0.081, 0.469, 0.488, 0.287, 0.14, 0.481, 0.769, 0.765, 0.207, 0.483, 0.559, 0.659, 0.147, 0.523, 0.136, 0.753, 0.509, 0.413, 0.546, 0.21, 0.866]
global q = [0.241, 0.741, 0.997, 0.644, 0.539, 0.975, 0.634, 0.38, 0.602, 0.533, 0.954, 0.916, 0.774, 0.943, 0.307, 0.935, 0.592, 0.989, 0.383, 0.404, 0.897, 0.958, 0.469, 0.962, 0.439, 0.721, 0.91, 0.848, 0.774, 0.768, 0.156, 0.899, 0.478, 0.648, 0.886, 0.479, 0.699, 0.306, 0.608, 0.998, 0.969, 0.943, 0.598, 0.915, 0.998, 0.481, 0.677, 0.904, 0.498, 0.963, 0.993, 0.394, 0.837, 0.895, 0.908, 0.767, 0.899, 0.877, 0.353, 0.775, 0.798, 0.874, 0.861, 0.664, 0.496, 0.997, 0.959, 0.92, 0.951, 0.913, 0.77, 0.951, 0.936, 0.898, 0.56, 0.645, 0.999, 0.377, 0.963, 0.791, 0.922, 0.829, 0.842, 0.813, 0.681, 0.837, 0.549, 0.954, 0.878, 0.885, 0.54, 0.649, 0.948, 0.862, 0.52, 0.981, 0.74, 0.838, 0.86, 0.609, 0.968, 0.945, 0.745, 0.973, 0.425, 0.899, 0.662, 0.769, 0.829, 0.15, 0.571, 0.94, 0.928, 0.787, 0.519, 0.902, 0.902, 0.998, 0.725, 0.824, 0.723, 0.759, 0.934, 0.999, 0.868, 0.879, 0.999, 0.957, 0.374, 0.633, 0.886, 0.391, 0.904, 0.806, 0.942, 0.806, 0.522, 0.969, 0.745, 0.715, 0.994, 0.826, 0.999, 0.479, 0.953, 0.401, 0.886, 0.218, 0.708, 0.88, 0.678, 0.954, 0.941, 0.852, 0.92, 0.349, 0.827, 0.893, 0.877, 0.883, 0.974, 0.692, 0.898, 0.702, 0.949, 0.716, 0.739, 0.882, 0.696, 0.907, 0.946, 0.979, 0.725, 0.858, 0.945, 0.279, 0.842, 0.76, 0.943, 0.88, 0.958, 0.656, 0.577, 0.803, 0.178, 0.519, 0.812, 0.453, 0.382, 0.274, 0.614, 0.884, 0.949, 0.828, 0.887, 0.983, 0.948, 0.989, 0.725, 0.789, 0.815, 0.892, 0.937, 0.786, 0.969, 0.02, 0.19, 0.528, 0.987, 0.651, 0.8, 0.989, 0.471, 0.833, 0.59, 0.932, 0.701, 0.533, 0.877, 0.854, 0.251, 0.927, 0.972, 0.942, 0.951, 0.926, 0.881, 0.921, 0.842, 0.712, 0.966, 0.907, 0.982, 0.964, 0.625, 0.593, 0.488, 0.39, 0.765, 0.778, 0.971, 0.872, 0.659, 0.804, 0.927, 0.234, 0.538, 0.672, 0.725, 0.964, 0.628, 0.947, 0.773, 0.773, 0.927, 0.923, 0.768, 0.209, 0.724, 0.442, 0.756, 0.719, 0.767, 0.844, 0.383, 0.879]
global origin = 1
global destination = 50