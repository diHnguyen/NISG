global arcs = [1 8; 1 20; 2 17; 2 23; 2 26; 3 9; 3 14; 3 20; 3 31; 3 35; 3 37; 4 2; 4 13; 4 19; 4 22; 4 26; 4 28; 4 37; 4 40; 5 6; 5 12; 5 17; 5 24; 5 30; 5 33; 6 9; 6 10; 6 30; 6 35; 7 10; 7 24; 7 33; 8 6; 8 14; 8 22; 8 31; 8 39; 9 4; 9 14; 9 29; 9 34; 10 5; 10 20; 10 22; 10 38; 11 10; 11 21; 11 40; 12 2; 12 10; 12 15; 12 18; 12 20; 12 24; 12 33; 12 34; 13 9; 13 26; 14 26; 14 29; 14 33; 14 35; 14 36; 15 9; 15 20; 15 25; 15 30; 15 31; 16 12; 16 19; 16 21; 16 25; 17 18; 17 19; 17 36; 17 39; 18 23; 18 24; 18 28; 18 40; 19 29; 20 28; 21 12; 21 15; 21 23; 21 29; 21 39; 22 3; 22 11; 22 13; 22 23; 23 4; 23 14; 23 30; 24 20; 25 6; 25 11; 25 18; 25 21; 25 29; 25 35; 26 6; 26 9; 26 11; 26 24; 27 25; 27 30; 27 32; 27 33; 28 2; 28 17; 28 21; 28 34; 28 35; 29 4; 29 13; 29 26; 30 10; 30 16; 30 19; 30 21; 30 31; 31 6; 31 8; 31 27; 32 20; 32 22; 32 23; 32 25; 32 31; 32 38; 33 4; 33 7; 33 23; 33 24; 33 39; 34 3; 34 17; 34 18; 35 9; 35 19; 35 26; 36 14; 36 22; 36 24; 36 29; 36 33; 36 40; 37 4; 37 10; 37 17; 37 19; 37 27; 37 29; 37 38; 38 26; 38 29; 38 32; 38 35; 39 6; 39 7; 39 9]
global d_x = [2.0, 4.0, 7.0, 6.0, 7.0, 8.0, 9.0, 7.0, 6.0, 2.0, 4.0, 7.0, 1.0, 9.0, 4.0, 6.0, 5.0, 3.0, 1.0, 4.0, 4.0, 9.0, 9.0, 6.0, 4.0, 4.0, 10.0, 9.0, 4.0, 5.0, 9.0, 9.0, 3.0, 8.0, 6.0, 5.0, 6.0, 7.0, 3.0, 1.0, 2.0, 10.0, 2.0, 4.0, 3.0, 10.0, 5.0, 3.0, 7.0, 8.0, 2.0, 3.0, 2.0, 3.0, 4.0, 4.0, 6.0, 6.0, 3.0, 5.0, 1.0, 10.0, 1.0, 3.0, 4.0, 8.0, 2.0, 9.0, 8.0, 4.0, 9.0, 5.0, 5.0, 7.0, 4.0, 2.0, 10.0, 2.0, 8.0, 5.0, 9.0, 4.0, 6.0, 7.0, 7.0, 9.0, 4.0, 10.0, 8.0, 2.0, 9.0, 3.0, 7.0, 5.0, 2.0, 7.0, 7.0, 7.0, 6.0, 6.0, 8.0, 1.0, 5.0, 9.0, 10.0, 6.0, 7.0, 8.0, 6.0, 5.0, 2.0, 1.0, 8.0, 8.0, 2.0, 6.0, 4.0, 10.0, 7.0, 1.0, 1.0, 10.0, 3.0, 5.0, 4.0, 4.0, 2.0, 10.0, 3.0, 7.0, 6.0, 5.0, 1.0, 4.0, 9.0, 9.0, 7.0, 1.0, 7.0, 4.0, 2.0, 5.0, 1.0, 5.0, 5.0, 6.0, 4.0, 4.0, 5.0, 4.0, 8.0, 8.0, 1.0, 3.0, 10.0, 4.0, 8.0, 2.0, 2.0, 4.0, 10.0, 1.0]
global b_x = 5
global d_y = [10.0, 5.0, 3.0, 1.0, 6.0, 1.0, 7.0, 9.0, 5.0, 2.0, 7.0, 6.0, 1.0, 4.0, 8.0, 5.0, 7.0, 6.0, 3.0, 8.0, 8.0, 6.0, 4.0, 6.0, 9.0, 2.0, 7.0, 5.0, 3.0, 2.0, 1.0, 1.0, 2.0, 4.0, 8.0, 4.0, 9.0, 3.0, 6.0, 2.0, 8.0, 4.0, 5.0, 4.0, 6.0, 8.0, 8.0, 5.0, 3.0, 10.0, 1.0, 9.0, 7.0, 4.0, 10.0, 9.0, 9.0, 7.0, 1.0, 7.0, 5.0, 3.0, 7.0, 4.0, 6.0, 2.0, 8.0, 7.0, 6.0, 6.0, 2.0, 4.0, 5.0, 1.0, 2.0, 4.0, 10.0, 3.0, 9.0, 6.0, 1.0, 5.0, 1.0, 4.0, 4.0, 6.0, 9.0, 1.0, 6.0, 6.0, 1.0, 7.0, 3.0, 7.0, 5.0, 4.0, 5.0, 10.0, 9.0, 1.0, 1.0, 4.0, 4.0, 5.0, 6.0, 3.0, 2.0, 10.0, 3.0, 9.0, 5.0, 10.0, 5.0, 4.0, 6.0, 5.0, 7.0, 1.0, 9.0, 9.0, 1.0, 2.0, 4.0, 8.0, 4.0, 8.0, 1.0, 2.0, 6.0, 6.0, 8.0, 4.0, 6.0, 10.0, 8.0, 4.0, 9.0, 10.0, 2.0, 7.0, 4.0, 10.0, 6.0, 10.0, 4.0, 8.0, 5.0, 3.0, 9.0, 9.0, 8.0, 5.0, 7.0, 2.0, 5.0, 5.0, 10.0, 5.0, 8.0, 4.0, 7.0, 3.0]
global b_y = 10
global p = [0.482, 0.61, 0.965, 0.641, 0.187, 0.171, 0.901, 0.963, 0.098, 0.292, 0.516, 0.214, 0.033, 0.147, 0.265, 0.969, 0.072, 0.261, 0.531, 0.534, 0.509, 0.567, 0.731, 0.851, 0.063, 0.66, 0.956, 0.799, 0.615, 0.824, 0.682, 0.472, 0.284, 0.001, 0.103, 0.589, 0.929, 0.851, 0.47, 0.578, 0.396, 0.898, 0.664, 0.932, 0.383, 0.058, 0.505, 0.803, 0.119, 0.951, 0.482, 0.942, 0.58, 0.042, 0.272, 0.047, 0.077, 0.879, 0.834, 0.898, 0.58, 0.515, 0.968, 0.39, 0.686, 0.651, 0.419, 0.603, 0.228, 0.508, 0.102, 0.042, 0.785, 0.678, 0.791, 0.674, 0.517, 0.944, 0.342, 0.847, 0.584, 0.726, 0.889, 0.536, 0.573, 0.224, 0.088, 0.667, 0.967, 0.472, 0.962, 0.378, 0.039, 0.351, 0.066, 0.643, 0.711, 0.109, 0.475, 0.431, 0.585, 0.176, 0.925, 0.508, 0.92, 0.253, 0.833, 0.716, 0.425, 0.173, 0.413, 0.965, 0.829, 0.459, 0.158, 0.404, 0.569, 0.28, 0.119, 0.872, 0.502, 0.912, 0.346, 0.276, 0.537, 0.755, 0.13, 0.77, 0.832, 0.93, 0.185, 0.117, 0.552, 0.479, 0.382, 0.679, 0.555, 0.985, 0.309, 0.511, 0.595, 0.839, 0.639, 0.708, 0.08, 0.954, 0.14, 0.935, 0.833, 0.014, 0.979, 0.183, 0.613, 0.925, 0.966, 0.571, 0.074, 0.256, 0.476, 0.195, 0.84, 0.972]
global q = [0.643, 0.63, 0.987, 0.871, 0.694, 0.189, 0.949, 0.984, 0.195, 0.532, 0.606, 0.553, 0.728, 0.221, 0.575, 0.983, 0.198, 0.437, 0.855, 0.653, 0.932, 0.812, 0.812, 0.956, 0.286, 0.706, 0.972, 0.873, 0.813, 0.974, 0.893, 0.557, 0.988, 0.888, 0.891, 0.664, 0.979, 0.947, 0.656, 0.585, 0.65, 0.964, 0.914, 0.939, 0.497, 0.413, 0.662, 0.998, 0.883, 0.999, 0.492, 0.992, 0.746, 0.748, 0.77, 0.279, 0.735, 0.953, 0.948, 0.979, 0.714, 0.696, 0.988, 0.423, 0.897, 0.667, 0.502, 0.952, 0.348, 0.544, 0.971, 0.766, 0.802, 0.91, 0.819, 0.971, 0.994, 0.948, 0.644, 0.874, 0.725, 0.832, 0.916, 0.785, 0.608, 0.871, 0.851, 0.705, 0.991, 0.995, 0.992, 0.696, 0.059, 0.666, 0.182, 0.869, 0.888, 0.753, 0.511, 0.531, 0.951, 0.908, 0.961, 0.686, 0.986, 0.785, 0.847, 0.84, 0.679, 0.518, 0.804, 0.982, 0.942, 0.685, 0.84, 0.44, 0.805, 0.675, 0.489, 0.96, 0.73, 0.95, 0.426, 0.683, 0.638, 0.924, 0.628, 0.863, 0.892, 0.96, 0.188, 0.456, 0.601, 0.598, 0.464, 0.775, 0.987, 0.996, 0.833, 0.713, 0.976, 0.848, 0.851, 0.749, 0.358, 0.991, 0.439, 0.99, 0.875, 0.554, 0.993, 0.215, 0.639, 0.959, 0.969, 0.987, 0.232, 0.442, 0.975, 0.723, 0.845, 0.997]
global origin = 1
global destination = 40
