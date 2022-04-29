global arcs = [1 15; 1 18; 1 23; 1 32; 2 13; 2 17; 2 22; 3 11; 3 34; 4 28; 5 8; 5 9; 5 25; 5 26; 5 32; 5 34; 6 5; 6 15; 6 27; 6 33; 7 3; 7 5; 7 12; 7 14; 7 29; 7 35; 8 7; 8 15; 8 17; 9 3; 9 8; 9 11; 9 15; 9 30; 10 3; 10 14; 10 19; 10 20; 10 29; 11 15; 11 18; 11 26; 11 27; 11 35; 12 20; 12 24; 12 25; 13 15; 13 19; 13 22; 13 24; 14 16; 14 31; 14 33; 15 3; 15 12; 15 13; 15 19; 15 28; 15 34; 16 5; 16 22; 16 23; 16 34; 16 35; 17 3; 17 5; 17 6; 17 10; 17 11; 17 14; 17 16; 17 30; 17 31; 18 6; 18 13; 18 20; 18 34; 19 3; 19 22; 19 33; 20 2; 20 3; 20 8; 20 12; 20 14; 20 19; 20 23; 21 2; 21 3; 21 6; 21 28; 21 33; 22 12; 22 13; 22 20; 22 28; 23 13; 23 18; 23 30; 23 34; 24 6; 24 16; 24 22; 24 33; 25 3; 25 11; 25 26; 25 29; 25 34; 26 6; 26 23; 27 6; 27 19; 28 8; 28 33; 29 4; 29 6; 29 8; 29 14; 29 17; 29 26; 29 31; 30 3; 30 4; 30 6; 31 9; 31 15; 31 17; 31 34; 32 34; 33 21; 33 25; 34 3; 34 25]
global d_x = [2.0, 1.0, 7.0, 6.0, 4.0, 3.0, 9.0, 10.0, 6.0, 1.0, 9.0, 9.0, 1.0, 9.0, 4.0, 10.0, 4.0, 9.0, 7.0, 7.0, 1.0, 7.0, 10.0, 7.0, 8.0, 7.0, 2.0, 9.0, 10.0, 3.0, 6.0, 5.0, 1.0, 5.0, 5.0, 5.0, 9.0, 1.0, 2.0, 2.0, 3.0, 9.0, 9.0, 8.0, 6.0, 5.0, 1.0, 10.0, 5.0, 10.0, 3.0, 3.0, 4.0, 8.0, 4.0, 6.0, 2.0, 5.0, 10.0, 1.0, 6.0, 3.0, 5.0, 10.0, 6.0, 2.0, 1.0, 10.0, 9.0, 5.0, 7.0, 2.0, 9.0, 4.0, 4.0, 3.0, 4.0, 7.0, 8.0, 1.0, 7.0, 5.0, 7.0, 10.0, 2.0, 8.0, 9.0, 3.0, 8.0, 10.0, 8.0, 1.0, 3.0, 4.0, 8.0, 1.0, 5.0, 1.0, 10.0, 9.0, 6.0, 1.0, 8.0, 4.0, 8.0, 2.0, 1.0, 5.0, 3.0, 2.0, 7.0, 8.0, 6.0, 3.0, 6.0, 5.0, 8.0, 6.0, 10.0, 1.0, 2.0, 2.0, 4.0, 5.0, 2.0, 3.0, 1.0, 5.0, 6.0, 10.0, 3.0, 8.0, 8.0, 2.0, 4.0]
global b_x = 5
global d_y = [5.0, 5.0, 1.0, 8.0, 6.0, 5.0, 7.0, 3.0, 10.0, 9.0, 8.0, 10.0, 4.0, 10.0, 7.0, 7.0, 3.0, 7.0, 7.0, 3.0, 8.0, 6.0, 6.0, 9.0, 9.0, 10.0, 9.0, 5.0, 4.0, 8.0, 1.0, 10.0, 9.0, 6.0, 4.0, 6.0, 4.0, 6.0, 5.0, 2.0, 7.0, 2.0, 4.0, 3.0, 5.0, 1.0, 10.0, 1.0, 7.0, 6.0, 4.0, 5.0, 7.0, 1.0, 6.0, 9.0, 4.0, 4.0, 7.0, 9.0, 9.0, 7.0, 5.0, 10.0, 10.0, 5.0, 10.0, 10.0, 3.0, 5.0, 7.0, 6.0, 5.0, 3.0, 7.0, 9.0, 7.0, 3.0, 3.0, 5.0, 7.0, 6.0, 2.0, 5.0, 7.0, 4.0, 4.0, 6.0, 4.0, 8.0, 8.0, 8.0, 5.0, 5.0, 1.0, 8.0, 3.0, 7.0, 3.0, 7.0, 5.0, 2.0, 2.0, 2.0, 8.0, 6.0, 5.0, 4.0, 7.0, 4.0, 8.0, 4.0, 6.0, 4.0, 8.0, 10.0, 5.0, 3.0, 1.0, 8.0, 3.0, 8.0, 9.0, 2.0, 1.0, 10.0, 10.0, 6.0, 5.0, 2.0, 5.0, 2.0, 7.0, 1.0, 2.0]
global b_y = 10
global p = [0.329, 0.946, 0.215, 0.288, 0.298, 0.962, 0.035, 0.553, 0.233, 0.424, 0.94, 0.814, 0.009, 0.373, 0.384, 0.518, 0.974, 0.8, 0.732, 0.647, 0.299, 0.805, 0.631, 0.143, 0.252, 0.891, 0.509, 0.518, 0.71, 0.53, 0.377, 0.396, 0.383, 0.778, 0.617, 0.655, 0.377, 0.545, 0.645, 0.73, 0.542, 0.989, 0.949, 0.642, 0.173, 0.525, 0.095, 0.778, 0.658, 0.338, 0.706, 0.284, 0.38, 0.356, 0.579, 0.148, 0.764, 0.295, 0.244, 0.322, 0.907, 0.407, 0.092, 0.612, 0.32, 0.197, 0.441, 0.485, 0.263, 0.986, 0.883, 0.208, 0.227, 0.471, 0.213, 0.827, 0.999, 0.265, 0.338, 0.278, 0.406, 0.41, 0.167, 0.562, 0.869, 0.958, 0.824, 0.244, 0.351, 0.716, 0.55, 0.063, 0.187, 0.404, 0.728, 0.184, 0.648, 0.547, 0.363, 0.764, 0.646, 0.565, 0.228, 0.631, 0.071, 0.668, 0.384, 0.936, 0.526, 0.491, 0.105, 0.552, 0.532, 0.431, 0.512, 0.307, 0.081, 0.122, 0.012, 0.41, 0.703, 0.454, 0.533, 0.154, 0.305, 0.649, 0.299, 0.961, 0.914, 0.848, 0.465, 0.622, 0.769, 0.542, 0.672]
global q = [0.543, 0.993, 0.824, 0.493, 0.911, 0.972, 0.405, 0.729, 0.671, 0.424, 0.995, 0.859, 0.87, 0.637, 0.634, 0.694, 0.991, 0.982, 0.78, 0.722, 0.415, 0.935, 0.725, 0.98, 0.381, 0.995, 0.927, 0.557, 0.727, 0.539, 0.789, 0.659, 0.788, 0.859, 0.9, 0.679, 0.721, 0.629, 0.69, 0.772, 0.866, 0.991, 0.999, 0.724, 0.516, 0.569, 0.15, 0.809, 0.971, 0.497, 0.851, 0.627, 0.792, 0.414, 0.99, 0.643, 0.984, 0.583, 0.33, 0.889, 0.976, 0.941, 0.46, 0.844, 0.728, 0.706, 0.785, 0.836, 0.71, 0.987, 0.927, 0.841, 0.865, 0.826, 0.938, 0.933, 0.999, 0.793, 0.623, 0.617, 0.815, 0.833, 0.737, 0.591, 0.896, 0.986, 0.927, 0.688, 0.408, 0.833, 0.986, 0.666, 0.837, 0.469, 0.843, 0.755, 0.672, 0.555, 0.459, 0.836, 0.764, 0.643, 0.46, 0.768, 0.686, 0.918, 0.418, 0.939, 0.88, 0.662, 0.92, 0.812, 0.841, 0.529, 0.751, 0.794, 0.59, 0.747, 0.772, 0.891, 0.843, 0.853, 0.985, 0.972, 0.534, 0.841, 0.902, 0.986, 0.926, 0.848, 0.571, 0.726, 0.955, 0.934, 0.711]
global origin = 1
global destination = 35
