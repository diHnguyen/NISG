global arcs = [1 13; 1 18; 1 20; 1 27; 2 4; 2 9; 2 21; 2 22; 2 27; 3 4; 3 11; 3 23; 4 3; 4 7; 4 9; 4 17; 4 19; 5 19; 5 24; 5 39; 6 5; 6 17; 6 26; 7 28; 7 38; 8 39; 9 35; 10 12; 10 21; 10 23; 10 25; 10 27; 11 13; 11 17; 11 24; 11 30; 11 31; 11 33; 11 34; 11 36; 11 39; 12 2; 12 8; 12 18; 13 26; 13 30; 14 4; 14 7; 14 22; 14 29; 14 36; 14 38; 15 4; 15 19; 15 26; 15 35; 15 39; 16 21; 16 26; 16 32; 17 6; 17 13; 17 14; 17 22; 17 28; 18 16; 18 34; 18 35; 18 37; 19 7; 19 20; 19 26; 20 2; 20 37; 20 39; 21 12; 21 33; 22 3; 23 13; 23 24; 24 6; 24 13; 24 32; 24 33; 24 34; 24 40; 25 5; 25 19; 25 29; 25 40; 26 15; 26 20; 26 33; 27 15; 27 20; 27 40; 28 22; 28 27; 28 33; 28 35; 29 3; 29 12; 29 19; 29 24; 29 30; 29 40; 30 3; 30 31; 30 33; 31 4; 31 11; 31 20; 32 10; 32 36; 33 39; 34 17; 34 28; 34 29; 35 7; 35 9; 35 14; 35 27; 35 40; 36 3; 36 32; 37 2; 37 7; 37 12; 37 13; 37 14; 37 21; 37 33; 38 25; 38 27; 38 32; 39 11; 39 16; 39 31; 39 33; 39 35; 39 40]
global d_x = [5.0, 1.0, 1.0, 7.0, 6.0, 5.0, 1.0, 4.0, 6.0, 10.0, 2.0, 1.0, 4.0, 2.0, 7.0, 6.0, 2.0, 4.0, 6.0, 7.0, 2.0, 9.0, 1.0, 4.0, 6.0, 6.0, 1.0, 4.0, 2.0, 4.0, 5.0, 5.0, 3.0, 2.0, 9.0, 5.0, 1.0, 3.0, 6.0, 4.0, 3.0, 1.0, 5.0, 3.0, 8.0, 3.0, 3.0, 8.0, 7.0, 7.0, 4.0, 5.0, 1.0, 6.0, 9.0, 1.0, 9.0, 1.0, 1.0, 1.0, 6.0, 4.0, 8.0, 3.0, 5.0, 4.0, 7.0, 6.0, 7.0, 2.0, 1.0, 9.0, 6.0, 2.0, 1.0, 8.0, 3.0, 7.0, 2.0, 7.0, 5.0, 4.0, 2.0, 7.0, 3.0, 9.0, 6.0, 3.0, 10.0, 6.0, 5.0, 3.0, 2.0, 3.0, 10.0, 8.0, 7.0, 3.0, 1.0, 9.0, 6.0, 2.0, 8.0, 2.0, 5.0, 4.0, 4.0, 8.0, 10.0, 8.0, 2.0, 8.0, 2.0, 7.0, 8.0, 10.0, 4.0, 10.0, 10.0, 9.0, 10.0, 7.0, 2.0, 2.0, 7.0, 4.0, 5.0, 2.0, 2.0, 5.0, 2.0, 2.0, 7.0, 9.0, 4.0, 3.0, 5.0, 1.0, 5.0, 6.0, 4.0]
global b_x = 5
global d_y = [5.0, 9.0, 7.0, 2.0, 4.0, 7.0, 2.0, 7.0, 5.0, 10.0, 7.0, 1.0, 10.0, 4.0, 10.0, 8.0, 8.0, 2.0, 2.0, 5.0, 3.0, 7.0, 6.0, 10.0, 7.0, 10.0, 6.0, 9.0, 3.0, 6.0, 8.0, 6.0, 8.0, 1.0, 7.0, 10.0, 8.0, 8.0, 5.0, 8.0, 9.0, 9.0, 3.0, 10.0, 3.0, 3.0, 1.0, 2.0, 5.0, 9.0, 7.0, 6.0, 10.0, 5.0, 2.0, 1.0, 5.0, 3.0, 7.0, 3.0, 4.0, 1.0, 9.0, 2.0, 1.0, 8.0, 8.0, 8.0, 3.0, 5.0, 3.0, 8.0, 7.0, 9.0, 5.0, 8.0, 6.0, 5.0, 2.0, 10.0, 7.0, 5.0, 5.0, 5.0, 6.0, 4.0, 8.0, 1.0, 2.0, 1.0, 8.0, 7.0, 5.0, 10.0, 7.0, 8.0, 7.0, 5.0, 5.0, 6.0, 9.0, 6.0, 10.0, 3.0, 8.0, 9.0, 5.0, 8.0, 4.0, 10.0, 2.0, 5.0, 10.0, 10.0, 1.0, 4.0, 10.0, 1.0, 9.0, 4.0, 8.0, 10.0, 9.0, 5.0, 5.0, 9.0, 3.0, 5.0, 9.0, 7.0, 9.0, 5.0, 10.0, 6.0, 1.0, 10.0, 5.0, 6.0, 3.0, 2.0, 7.0]
global b_y = 10
global p = [0.242, 0.388, 0.987, 0.063, 0.169, 0.06, 0.078, 0.778, 0.265, 0.78, 0.181, 0.392, 0.824, 0.628, 0.804, 0.345, 0.144, 0.573, 0.424, 0.555, 0.39, 0.599, 0.092, 0.625, 0.723, 0.352, 0.534, 0.972, 0.475, 0.48, 0.746, 0.059, 0.051, 0.469, 0.677, 0.506, 0.885, 0.785, 0.473, 0.582, 0.73, 0.228, 0.205, 0.469, 0.229, 0.307, 0.009, 0.954, 0.385, 0.005, 0.826, 0.644, 0.743, 0.027, 0.271, 0.472, 0.826, 0.761, 0.069, 0.572, 0.283, 0.046, 0.922, 0.815, 0.363, 0.332, 0.007, 0.319, 0.122, 0.217, 0.018, 0.119, 0.254, 0.414, 0.645, 0.88, 0.681, 0.183, 0.69, 0.151, 0.789, 0.871, 0.9, 0.874, 0.249, 0.631, 0.263, 0.967, 0.502, 0.067, 0.166, 0.352, 0.288, 0.408, 0.639, 0.753, 0.766, 0.552, 0.787, 0.814, 0.335, 0.053, 0.97, 0.289, 0.752, 0.487, 0.489, 0.476, 0.987, 0.579, 0.404, 0.071, 0.842, 0.103, 0.23, 0.961, 0.234, 0.635, 0.848, 0.714, 0.041, 0.533, 0.284, 0.248, 0.511, 0.37, 0.817, 0.027, 0.726, 0.048, 0.965, 0.319, 0.139, 0.854, 0.957, 0.461, 0.911, 0.171, 0.647, 0.971, 0.264]
global q = [0.375, 0.933, 0.992, 0.532, 0.417, 0.907, 0.112, 0.834, 0.387, 0.828, 0.759, 0.992, 0.99, 0.883, 0.973, 0.5, 0.212, 0.841, 0.84, 0.673, 0.69, 0.998, 0.501, 0.986, 0.872, 0.763, 0.968, 0.998, 0.706, 0.592, 0.823, 0.186, 0.628, 0.663, 0.87, 0.586, 0.973, 0.91, 0.848, 0.698, 0.812, 0.692, 0.705, 0.487, 0.794, 0.87, 0.634, 0.999, 0.756, 0.186, 0.883, 0.739, 0.963, 0.826, 0.864, 0.87, 0.828, 0.831, 0.953, 0.812, 0.673, 0.058, 0.962, 0.957, 0.908, 0.886, 0.653, 0.757, 0.155, 0.253, 0.1, 0.181, 0.402, 0.917, 0.944, 0.956, 0.791, 0.76, 0.703, 0.241, 0.84, 0.929, 0.935, 0.88, 0.49, 0.872, 0.859, 0.972, 0.699, 0.45, 0.744, 0.459, 0.312, 0.842, 0.987, 0.966, 0.851, 0.792, 0.976, 0.941, 0.436, 0.908, 0.99, 0.761, 0.839, 0.607, 0.986, 0.585, 0.993, 0.77, 0.744, 0.333, 0.919, 0.382, 0.41, 0.984, 0.896, 0.755, 0.926, 0.81, 0.091, 0.743, 0.821, 0.449, 0.969, 0.836, 0.922, 0.217, 0.861, 0.261, 0.984, 0.468, 0.396, 0.986, 0.957, 0.804, 0.918, 0.877, 0.884, 0.975, 0.6]
global origin = 1
global destination = 40