global arcs = [1 26; 1 35; 1 37; 2 6; 2 21; 2 32; 3 4; 3 12; 3 25; 3 39; 4 39; 5 3; 5 16; 5 23; 5 27; 6 3; 6 8; 6 18; 6 29; 6 35; 7 11; 7 12; 7 20; 7 21; 7 32; 7 33; 7 36; 8 5; 8 10; 8 33; 9 6; 9 7; 9 24; 9 37; 9 38; 10 8; 10 9; 10 11; 10 14; 10 17; 10 26; 10 35; 10 40; 11 2; 11 7; 11 13; 11 15; 11 31; 12 7; 12 34; 12 39; 12 40; 13 20; 13 37; 14 9; 14 11; 15 14; 15 19; 16 7; 16 9; 16 28; 17 22; 17 36; 17 40; 18 14; 18 26; 18 35; 19 15; 19 22; 19 26; 19 35; 20 30; 20 33; 21 6; 21 25; 21 30; 21 40; 22 3; 22 19; 23 15; 23 27; 23 30; 23 32; 23 33; 24 3; 24 12; 24 33; 25 6; 25 14; 25 19; 25 20; 26 15; 26 25; 26 36; 27 3; 27 20; 27 23; 27 33; 27 34; 27 37; 27 38; 28 8; 28 9; 28 24; 28 33; 28 39; 29 22; 29 23; 29 26; 29 31; 29 36; 30 19; 30 21; 30 36; 31 18; 32 5; 32 6; 32 26; 32 27; 32 28; 33 12; 33 24; 34 6; 34 10; 34 29; 34 35; 34 37; 34 38; 35 18; 35 37; 35 38; 36 4; 36 8; 36 15; 36 16; 36 21; 37 3; 37 10; 37 11; 37 18; 37 38; 38 15; 38 17; 39 15]
global d_x = [5.0, 10.0, 4.0, 8.0, 2.0, 4.0, 10.0, 1.0, 6.0, 8.0, 4.0, 9.0, 8.0, 8.0, 5.0, 3.0, 9.0, 2.0, 5.0, 10.0, 5.0, 4.0, 3.0, 1.0, 6.0, 6.0, 5.0, 4.0, 9.0, 6.0, 7.0, 6.0, 3.0, 10.0, 4.0, 3.0, 7.0, 5.0, 10.0, 10.0, 1.0, 7.0, 9.0, 10.0, 8.0, 9.0, 6.0, 7.0, 2.0, 2.0, 8.0, 10.0, 8.0, 7.0, 6.0, 10.0, 3.0, 10.0, 5.0, 7.0, 10.0, 1.0, 1.0, 4.0, 7.0, 3.0, 9.0, 6.0, 1.0, 2.0, 3.0, 10.0, 4.0, 4.0, 8.0, 2.0, 1.0, 2.0, 7.0, 4.0, 9.0, 3.0, 10.0, 1.0, 1.0, 3.0, 1.0, 4.0, 10.0, 2.0, 7.0, 10.0, 2.0, 3.0, 7.0, 8.0, 6.0, 5.0, 10.0, 7.0, 3.0, 9.0, 2.0, 1.0, 1.0, 5.0, 7.0, 6.0, 5.0, 3.0, 2.0, 3.0, 6.0, 2.0, 1.0, 4.0, 3.0, 9.0, 3.0, 5.0, 7.0, 4.0, 5.0, 10.0, 4.0, 7.0, 2.0, 6.0, 1.0, 5.0, 5.0, 7.0, 6.0, 10.0, 3.0, 9.0, 6.0, 5.0, 7.0, 1.0, 5.0, 4.0, 10.0, 2.0]
global b_x = 5
global d_y = [8.0, 8.0, 10.0, 1.0, 9.0, 6.0, 1.0, 9.0, 9.0, 6.0, 3.0, 3.0, 7.0, 9.0, 1.0, 6.0, 2.0, 10.0, 6.0, 10.0, 2.0, 3.0, 7.0, 4.0, 5.0, 8.0, 7.0, 3.0, 7.0, 1.0, 7.0, 8.0, 5.0, 6.0, 9.0, 7.0, 10.0, 3.0, 7.0, 4.0, 4.0, 7.0, 5.0, 7.0, 7.0, 6.0, 6.0, 2.0, 1.0, 7.0, 9.0, 1.0, 4.0, 6.0, 3.0, 7.0, 4.0, 8.0, 9.0, 6.0, 4.0, 1.0, 1.0, 3.0, 7.0, 10.0, 2.0, 8.0, 8.0, 7.0, 10.0, 9.0, 10.0, 4.0, 1.0, 2.0, 5.0, 5.0, 2.0, 8.0, 5.0, 2.0, 1.0, 2.0, 3.0, 1.0, 6.0, 7.0, 5.0, 9.0, 2.0, 8.0, 6.0, 10.0, 7.0, 7.0, 5.0, 1.0, 7.0, 2.0, 5.0, 1.0, 4.0, 2.0, 10.0, 6.0, 8.0, 2.0, 1.0, 5.0, 8.0, 6.0, 10.0, 3.0, 8.0, 6.0, 1.0, 6.0, 9.0, 5.0, 3.0, 8.0, 9.0, 3.0, 4.0, 2.0, 2.0, 4.0, 6.0, 7.0, 3.0, 5.0, 10.0, 8.0, 3.0, 5.0, 9.0, 1.0, 9.0, 6.0, 9.0, 8.0, 1.0, 2.0]
global b_y = 10
global p = [0.879, 0.65, 0.379, 0.87, 0.076, 0.233, 0.662, 0.57, 0.535, 0.461, 0.977, 0.252, 0.149, 0.961, 0.924, 0.967, 0.887, 0.654, 0.486, 0.234, 0.376, 0.88, 0.865, 0.348, 0.28, 0.018, 0.481, 0.239, 0.316, 0.304, 0.324, 0.709, 0.214, 0.859, 0.19, 0.243, 0.758, 0.118, 0.437, 0.858, 0.471, 0.535, 0.571, 0.86, 0.885, 0.943, 0.055, 0.262, 0.744, 0.879, 0.301, 0.157, 0.83, 0.85, 0.696, 0.415, 0.89, 0.773, 0.674, 0.849, 0.492, 0.745, 0.054, 0.72, 0.863, 0.342, 0.615, 0.622, 0.513, 0.528, 0.375, 0.16, 0.472, 0.66, 0.403, 0.399, 0.879, 0.471, 0.572, 0.036, 0.428, 0.88, 0.746, 0.82, 0.489, 0.454, 0.628, 0.228, 0.292, 0.195, 0.839, 0.93, 0.719, 0.165, 0.681, 0.145, 0.042, 0.028, 0.537, 0.151, 0.668, 0.021, 0.769, 0.398, 0.761, 0.185, 0.064, 0.983, 0.119, 0.055, 0.224, 0.73, 0.998, 0.194, 0.881, 0.08, 0.975, 0.13, 0.966, 0.898, 0.308, 0.795, 0.723, 0.624, 0.862, 0.293, 0.862, 0.533, 0.602, 0.867, 0.818, 0.916, 0.773, 0.163, 0.645, 0.038, 0.711, 0.724, 0.116, 0.628, 0.888, 0.712, 0.126, 0.889]
global q = [0.95, 0.696, 0.952, 0.96, 0.572, 0.292, 0.901, 0.872, 0.792, 0.681, 0.992, 0.751, 0.17, 0.991, 0.952, 0.969, 0.944, 0.844, 0.843, 0.93, 0.453, 0.901, 0.933, 0.649, 0.839, 0.996, 0.843, 0.486, 0.912, 0.957, 0.923, 0.81, 0.214, 0.859, 0.865, 0.546, 0.894, 0.316, 0.994, 0.99, 0.655, 0.612, 0.8, 0.927, 0.916, 0.943, 0.59, 0.585, 0.755, 0.957, 0.927, 0.67, 0.928, 0.88, 0.975, 0.606, 0.973, 0.99, 0.759, 0.946, 0.728, 0.875, 0.331, 0.952, 0.939, 0.932, 0.819, 0.703, 0.921, 0.627, 0.67, 0.943, 0.617, 0.825, 0.75, 0.661, 0.952, 0.634, 0.751, 0.151, 0.82, 0.988, 0.871, 0.846, 0.78, 0.46, 0.779, 0.877, 0.942, 0.502, 0.891, 0.969, 0.734, 0.466, 0.884, 0.92, 0.387, 0.788, 0.612, 0.632, 0.951, 0.773, 0.965, 0.589, 0.783, 0.38, 0.44, 0.993, 0.782, 0.964, 0.739, 0.952, 0.999, 0.691, 0.896, 0.689, 0.976, 0.975, 0.98, 0.898, 0.639, 0.835, 0.984, 0.706, 0.952, 0.535, 0.911, 0.732, 0.749, 0.894, 0.851, 0.952, 0.994, 0.986, 0.836, 0.809, 0.872, 0.963, 0.595, 0.663, 0.919, 0.778, 0.535, 0.924]
global origin = 1
global destination = 40
