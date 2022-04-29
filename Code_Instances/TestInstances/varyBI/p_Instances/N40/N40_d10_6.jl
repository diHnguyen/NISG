global arcs = [1 25; 2 16; 2 20; 2 22; 3 7; 3 12; 3 17; 3 39; 4 21; 4 22; 4 29; 4 34; 4 35; 4 38; 5 10; 5 11; 5 26; 5 30; 6 15; 6 21; 7 12; 7 33; 7 36; 8 5; 8 10; 8 24; 8 39; 9 12; 9 19; 9 20; 9 37; 9 38; 10 3; 10 7; 10 9; 10 17; 10 27; 10 28; 10 31; 10 34; 10 39; 11 5; 11 14; 11 20; 11 27; 11 33; 11 36; 12 28; 13 4; 13 6; 13 14; 13 25; 13 28; 13 34; 13 38; 14 9; 14 23; 15 26; 15 34; 16 12; 16 35; 16 37; 17 19; 17 22; 17 31; 17 39; 18 17; 18 29; 19 5; 19 27; 19 30; 19 33; 19 36; 19 38; 20 19; 20 27; 20 37; 21 17; 21 29; 22 2; 22 11; 22 25; 22 26; 22 33; 23 10; 23 14; 23 25; 23 27; 23 31; 24 3; 24 9; 24 19; 24 20; 24 34; 24 36; 25 7; 25 18; 26 8; 26 12; 26 38; 27 4; 27 29; 27 32; 28 4; 28 6; 28 31; 28 40; 29 11; 29 25; 29 30; 29 31; 29 39; 29 40; 30 4; 30 6; 30 28; 30 29; 30 31; 30 33; 31 20; 31 22; 31 34; 32 27; 32 34; 33 16; 33 18; 33 31; 34 3; 34 7; 34 21; 34 39; 35 37; 36 26; 36 29; 36 31; 37 2; 37 7; 37 10; 37 13; 37 21; 37 26; 37 30; 37 35; 38 11; 38 19; 38 22; 38 36; 39 4; 39 6; 39 11]
global d_x = [6.0, 8.0, 10.0, 1.0, 2.0, 1.0, 8.0, 10.0, 2.0, 10.0, 1.0, 10.0, 7.0, 8.0, 8.0, 3.0, 5.0, 1.0, 2.0, 5.0, 5.0, 2.0, 9.0, 2.0, 4.0, 6.0, 3.0, 10.0, 5.0, 7.0, 8.0, 7.0, 8.0, 9.0, 4.0, 9.0, 1.0, 3.0, 6.0, 6.0, 2.0, 8.0, 4.0, 10.0, 7.0, 2.0, 8.0, 7.0, 4.0, 7.0, 8.0, 9.0, 7.0, 5.0, 7.0, 10.0, 2.0, 6.0, 1.0, 4.0, 8.0, 3.0, 3.0, 7.0, 7.0, 10.0, 5.0, 9.0, 4.0, 6.0, 1.0, 3.0, 3.0, 3.0, 4.0, 6.0, 10.0, 9.0, 2.0, 4.0, 6.0, 1.0, 10.0, 2.0, 5.0, 2.0, 10.0, 6.0, 2.0, 9.0, 1.0, 4.0, 8.0, 8.0, 3.0, 7.0, 7.0, 2.0, 2.0, 7.0, 10.0, 1.0, 10.0, 2.0, 4.0, 9.0, 3.0, 5.0, 2.0, 2.0, 5.0, 2.0, 4.0, 1.0, 4.0, 7.0, 6.0, 3.0, 2.0, 7.0, 10.0, 6.0, 7.0, 2.0, 3.0, 1.0, 7.0, 6.0, 5.0, 10.0, 4.0, 6.0, 10.0, 3.0, 1.0, 10.0, 3.0, 7.0, 6.0, 3.0, 2.0, 1.0, 2.0, 9.0, 5.0, 8.0, 6.0, 5.0, 6.0, 2.0]
global b_x = 5
global d_y = [6.0, 4.0, 9.0, 3.0, 1.0, 9.0, 7.0, 7.0, 6.0, 2.0, 8.0, 7.0, 3.0, 1.0, 1.0, 3.0, 1.0, 6.0, 1.0, 4.0, 6.0, 9.0, 5.0, 9.0, 9.0, 10.0, 5.0, 8.0, 9.0, 8.0, 5.0, 6.0, 6.0, 10.0, 9.0, 1.0, 5.0, 3.0, 10.0, 9.0, 5.0, 1.0, 5.0, 6.0, 8.0, 7.0, 8.0, 2.0, 2.0, 5.0, 10.0, 7.0, 8.0, 3.0, 5.0, 2.0, 3.0, 5.0, 2.0, 5.0, 8.0, 1.0, 4.0, 7.0, 4.0, 8.0, 6.0, 2.0, 8.0, 9.0, 6.0, 9.0, 8.0, 8.0, 10.0, 4.0, 6.0, 4.0, 4.0, 8.0, 1.0, 7.0, 5.0, 6.0, 2.0, 2.0, 3.0, 6.0, 6.0, 10.0, 5.0, 4.0, 3.0, 3.0, 4.0, 8.0, 4.0, 9.0, 10.0, 4.0, 6.0, 2.0, 5.0, 4.0, 2.0, 6.0, 4.0, 10.0, 2.0, 10.0, 1.0, 2.0, 6.0, 6.0, 6.0, 3.0, 1.0, 5.0, 10.0, 3.0, 4.0, 3.0, 10.0, 3.0, 3.0, 9.0, 3.0, 2.0, 4.0, 4.0, 9.0, 8.0, 10.0, 5.0, 3.0, 4.0, 2.0, 6.0, 3.0, 10.0, 5.0, 3.0, 10.0, 2.0, 2.0, 1.0, 10.0, 3.0, 6.0, 3.0]
global b_y = 10
global p = [0.047, 0.774, 0.032, 0.239, 0.365, 0.084, 0.775, 0.015, 0.438, 0.516, 0.995, 0.912, 0.235, 0.637, 0.799, 0.61, 0.014, 0.828, 0.636, 0.73, 0.224, 0.004, 0.057, 0.472, 0.175, 0.363, 0.92, 0.164, 0.063, 0.102, 0.207, 0.639, 0.263, 0.578, 0.951, 0.393, 0.385, 0.546, 0.238, 0.521, 0.722, 0.249, 0.058, 0.808, 0.188, 0.086, 0.514, 0.812, 0.635, 0.397, 0.341, 0.993, 0.274, 0.945, 0.492, 0.917, 0.098, 0.702, 0.201, 0.017, 0.389, 0.891, 0.62, 0.14, 0.383, 0.976, 0.302, 0.986, 0.045, 0.16, 0.663, 0.808, 0.213, 0.765, 0.779, 0.27, 0.499, 0.312, 0.115, 0.497, 0.574, 0.188, 0.67, 0.82, 0.662, 0.583, 0.404, 0.338, 0.491, 0.931, 0.247, 0.765, 0.014, 0.065, 0.933, 0.015, 0.486, 0.512, 0.293, 0.609, 0.292, 0.713, 0.48, 0.891, 0.291, 0.361, 0.184, 0.638, 0.372, 0.506, 0.459, 0.588, 0.795, 0.806, 0.219, 0.097, 0.155, 0.609, 0.862, 0.109, 0.514, 0.832, 0.169, 0.09, 0.804, 0.423, 0.641, 0.696, 0.715, 0.965, 0.7, 0.269, 0.851, 0.733, 0.826, 0.719, 0.142, 0.149, 0.067, 0.253, 0.889, 0.51, 0.686, 0.915, 0.582, 0.327, 0.603, 0.989, 0.629, 0.33]
global q = [0.339, 0.833, 0.148, 0.37, 0.836, 0.942, 0.957, 0.1, 0.661, 0.877, 0.998, 0.983, 0.696, 0.776, 0.897, 0.825, 0.89, 0.913, 0.889, 0.986, 0.618, 0.569, 0.119, 0.864, 0.984, 0.831, 0.981, 0.412, 0.321, 0.338, 0.579, 0.668, 0.661, 0.747, 0.956, 0.512, 0.552, 0.686, 0.461, 0.775, 0.794, 0.275, 0.597, 0.969, 0.603, 0.166, 0.927, 0.982, 0.741, 0.975, 0.577, 0.997, 0.872, 0.978, 0.676, 0.931, 0.842, 0.887, 0.612, 0.063, 0.746, 0.894, 0.841, 0.766, 0.965, 0.994, 0.906, 0.996, 0.229, 0.945, 0.874, 0.888, 0.597, 0.785, 0.971, 0.719, 0.743, 0.859, 0.677, 0.948, 0.992, 0.838, 0.914, 0.94, 0.983, 0.716, 0.404, 0.776, 0.841, 0.964, 0.413, 0.925, 0.581, 0.388, 0.957, 0.623, 0.925, 0.754, 0.893, 0.81, 0.55, 0.777, 0.716, 0.971, 0.943, 0.885, 0.277, 0.727, 0.936, 0.508, 0.983, 0.917, 0.938, 0.897, 0.985, 0.866, 0.508, 0.844, 0.949, 0.392, 0.84, 0.959, 0.73, 0.236, 0.921, 0.706, 0.873, 0.777, 0.956, 0.979, 0.736, 0.331, 0.852, 0.927, 0.854, 0.996, 0.402, 0.762, 0.867, 0.633, 0.897, 0.526, 0.875, 0.978, 0.804, 0.859, 0.798, 0.998, 0.824, 0.605]
global origin = 1
global destination = 40
