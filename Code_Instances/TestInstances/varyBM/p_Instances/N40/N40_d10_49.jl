global arcs = [1 5; 1 7; 1 25; 1 29; 2 14; 2 23; 3 11; 3 16; 3 20; 3 23; 3 24; 4 3; 4 7; 4 13; 4 18; 4 24; 4 26; 4 30; 4 39; 5 9; 5 34; 5 40; 6 7; 6 13; 6 28; 6 29; 6 30; 6 31; 7 32; 7 36; 8 9; 8 18; 8 21; 8 24; 8 27; 8 30; 8 39; 9 24; 9 36; 10 8; 10 34; 11 9; 11 34; 12 6; 12 22; 12 32; 13 6; 13 30; 13 35; 14 3; 14 5; 14 21; 14 26; 14 39; 15 8; 15 11; 15 13; 15 18; 16 18; 16 31; 17 31; 17 35; 17 39; 18 11; 18 20; 18 23; 18 27; 18 29; 18 38; 18 40; 19 34; 19 35; 20 5; 20 21; 21 15; 21 25; 21 30; 21 38; 22 3; 22 17; 22 29; 22 37; 23 31; 24 3; 24 14; 24 30; 25 22; 26 21; 26 29; 27 7; 27 9; 27 28; 27 36; 27 37; 28 2; 28 5; 28 8; 28 11; 28 13; 28 17; 28 35; 29 8; 29 22; 29 23; 29 36; 30 12; 30 20; 30 33; 31 11; 31 15; 31 21; 31 33; 31 34; 32 14; 32 16; 33 19; 33 20; 33 35; 33 40; 34 3; 34 4; 35 11; 35 15; 35 19; 35 28; 36 14; 36 17; 36 18; 36 29; 36 30; 37 12; 37 14; 37 21; 37 24; 37 38; 37 39; 38 5; 38 11; 38 28; 39 2; 39 6; 39 10; 39 11; 39 23; 39 38]
global d_x = [3.0, 2.0, 4.0, 8.0, 1.0, 5.0, 10.0, 10.0, 10.0, 9.0, 1.0, 4.0, 5.0, 3.0, 7.0, 6.0, 1.0, 3.0, 10.0, 1.0, 4.0, 2.0, 4.0, 1.0, 3.0, 6.0, 10.0, 4.0, 7.0, 6.0, 2.0, 4.0, 1.0, 8.0, 2.0, 8.0, 6.0, 7.0, 8.0, 9.0, 9.0, 5.0, 2.0, 3.0, 4.0, 1.0, 7.0, 2.0, 1.0, 7.0, 5.0, 10.0, 10.0, 9.0, 8.0, 10.0, 7.0, 5.0, 6.0, 5.0, 7.0, 10.0, 3.0, 3.0, 1.0, 6.0, 5.0, 5.0, 6.0, 7.0, 10.0, 4.0, 1.0, 7.0, 9.0, 10.0, 5.0, 6.0, 2.0, 2.0, 5.0, 5.0, 6.0, 9.0, 1.0, 10.0, 6.0, 8.0, 3.0, 8.0, 3.0, 8.0, 4.0, 6.0, 10.0, 4.0, 9.0, 8.0, 8.0, 1.0, 7.0, 9.0, 5.0, 4.0, 5.0, 2.0, 8.0, 5.0, 1.0, 2.0, 9.0, 1.0, 9.0, 6.0, 3.0, 10.0, 8.0, 4.0, 9.0, 5.0, 1.0, 9.0, 5.0, 6.0, 7.0, 5.0, 10.0, 1.0, 4.0, 2.0, 3.0, 1.0, 3.0, 1.0, 5.0, 3.0, 5.0, 1.0, 10.0, 10.0, 9.0, 3.0, 8.0, 1.0, 3.0]
global b_x = 5
global d_y = [4.0, 5.0, 5.0, 9.0, 5.0, 7.0, 3.0, 1.0, 9.0, 7.0, 10.0, 5.0, 6.0, 1.0, 2.0, 2.0, 7.0, 3.0, 4.0, 6.0, 2.0, 10.0, 9.0, 1.0, 10.0, 7.0, 5.0, 5.0, 10.0, 2.0, 9.0, 9.0, 8.0, 1.0, 1.0, 2.0, 10.0, 4.0, 4.0, 3.0, 5.0, 8.0, 7.0, 1.0, 8.0, 1.0, 5.0, 7.0, 3.0, 4.0, 6.0, 4.0, 3.0, 2.0, 2.0, 8.0, 8.0, 9.0, 5.0, 9.0, 3.0, 5.0, 9.0, 2.0, 8.0, 6.0, 8.0, 8.0, 7.0, 3.0, 6.0, 3.0, 9.0, 9.0, 4.0, 1.0, 2.0, 8.0, 6.0, 7.0, 7.0, 7.0, 5.0, 1.0, 9.0, 5.0, 3.0, 7.0, 1.0, 4.0, 4.0, 5.0, 1.0, 1.0, 10.0, 9.0, 6.0, 6.0, 5.0, 4.0, 9.0, 8.0, 5.0, 6.0, 7.0, 6.0, 5.0, 2.0, 8.0, 3.0, 9.0, 5.0, 6.0, 6.0, 7.0, 6.0, 8.0, 10.0, 9.0, 2.0, 9.0, 5.0, 2.0, 3.0, 4.0, 8.0, 6.0, 3.0, 7.0, 7.0, 1.0, 6.0, 6.0, 10.0, 3.0, 8.0, 5.0, 7.0, 2.0, 8.0, 2.0, 7.0, 1.0, 7.0, 4.0]
global b_y = 10
global p = [0.079, 0.628, 0.054, 0.652, 0.045, 0.397, 0.895, 0.107, 0.898, 0.472, 0.012, 0.464, 0.256, 0.145, 0.122, 0.937, 0.359, 0.68, 0.826, 0.385, 0.574, 0.82, 0.562, 0.013, 0.558, 0.979, 0.221, 0.165, 0.789, 0.734, 0.772, 0.625, 0.157, 0.733, 0.483, 0.805, 0.444, 0.731, 0.62, 0.34, 0.86, 0.765, 0.757, 0.27, 0.239, 0.732, 0.508, 0.501, 0.537, 0.948, 0.871, 0.835, 0.555, 0.176, 0.898, 0.324, 0.579, 0.424, 0.093, 0.22, 0.391, 0.281, 0.783, 0.943, 0.892, 0.153, 0.368, 0.587, 0.876, 0.877, 0.332, 0.475, 0.529, 0.544, 0.426, 0.882, 0.99, 0.273, 0.981, 0.764, 0.907, 0.149, 0.633, 0.577, 0.351, 0.328, 0.997, 0.807, 0.465, 0.976, 0.711, 0.467, 0.426, 0.46, 0.475, 0.057, 0.965, 0.764, 0.092, 0.057, 0.745, 0.607, 0.647, 0.972, 0.118, 0.738, 0.014, 0.047, 0.958, 0.885, 0.155, 0.08, 0.758, 0.936, 0.531, 0.979, 0.492, 0.713, 0.01, 0.217, 0.818, 0.049, 0.354, 0.933, 0.049, 0.144, 0.397, 0.11, 0.178, 0.568, 0.045, 0.432, 0.675, 0.842, 0.702, 0.882, 0.061, 0.614, 0.101, 0.089, 0.272, 0.268, 0.033, 0.54, 0.605]
global q = [0.669, 0.78, 0.645, 0.799, 0.084, 0.698, 0.942, 0.623, 0.911, 0.955, 0.466, 0.648, 0.506, 0.603, 0.412, 0.97, 0.929, 0.803, 0.851, 0.907, 0.73, 0.824, 0.923, 0.39, 0.704, 0.992, 0.801, 0.575, 0.805, 0.878, 0.961, 0.957, 0.36, 0.822, 0.515, 0.833, 0.88, 0.79, 0.895, 0.634, 0.894, 0.945, 0.965, 0.922, 0.922, 0.756, 0.985, 0.557, 0.966, 0.968, 0.948, 0.951, 0.95, 0.394, 0.949, 0.669, 0.932, 0.665, 0.203, 0.759, 0.635, 0.779, 0.944, 0.993, 0.955, 0.968, 0.644, 0.836, 0.904, 0.908, 0.877, 0.617, 0.628, 0.852, 0.672, 0.934, 0.997, 0.809, 0.991, 0.933, 0.991, 0.729, 0.898, 0.817, 0.537, 0.436, 0.999, 0.885, 0.467, 0.979, 0.793, 0.571, 0.675, 0.632, 0.707, 0.373, 0.97, 0.84, 0.677, 0.338, 0.791, 0.949, 0.729, 0.991, 0.565, 0.96, 0.654, 0.118, 0.996, 0.958, 0.593, 0.352, 0.902, 0.961, 0.862, 0.998, 0.527, 0.882, 0.131, 0.349, 0.818, 0.331, 0.846, 0.946, 0.899, 0.653, 0.928, 0.553, 0.472, 0.595, 0.745, 0.84, 0.799, 0.851, 0.725, 0.936, 0.866, 0.743, 0.544, 0.41, 0.967, 0.896, 0.08, 0.966, 0.709]
global origin = 1
global destination = 40
