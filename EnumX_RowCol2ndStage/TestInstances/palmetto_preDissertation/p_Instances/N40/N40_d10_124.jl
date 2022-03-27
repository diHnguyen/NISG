global arcs = [1 15; 1 16; 1 23; 1 32; 1 39; 2 12; 2 17; 2 19; 2 28; 2 37; 3 4; 3 26; 3 29; 3 40; 4 17; 4 19; 4 27; 4 32; 4 34; 4 39; 5 7; 5 22; 5 25; 5 26; 5 31; 5 34; 6 4; 6 14; 6 20; 6 24; 6 25; 6 32; 7 3; 7 4; 7 9; 7 26; 8 17; 8 20; 8 26; 9 6; 9 10; 9 16; 10 12; 10 14; 10 27; 10 30; 10 38; 11 3; 11 4; 11 15; 11 24; 12 8; 12 11; 12 31; 13 12; 13 18; 13 22; 14 11; 14 15; 14 30; 14 39; 15 3; 15 10; 15 32; 15 40; 16 13; 16 15; 16 24; 16 31; 16 40; 17 11; 17 14; 17 24; 17 37; 17 39; 18 4; 18 10; 18 12; 18 25; 18 28; 19 20; 19 23; 19 37; 19 39; 20 8; 20 26; 20 36; 21 5; 21 8; 21 11; 21 13; 21 37; 21 38; 22 11; 22 18; 22 29; 23 3; 23 10; 23 12; 23 19; 23 34; 24 3; 24 10; 24 11; 24 13; 24 28; 24 35; 25 2; 25 4; 25 5; 25 6; 25 13; 25 32; 25 39; 25 40; 26 18; 26 23; 26 27; 26 30; 26 39; 27 11; 27 15; 27 26; 27 29; 27 33; 28 7; 28 9; 28 18; 28 30; 28 32; 29 6; 29 18; 29 40; 30 13; 30 24; 30 29; 31 17; 31 20; 31 28; 31 29; 31 39; 32 4; 32 14; 33 10; 33 14; 33 20; 33 21; 33 40; 34 21; 34 25; 35 8; 35 13; 35 17; 35 21; 35 36; 36 21; 36 37; 36 38; 37 6; 37 13; 37 17; 37 36; 37 39; 38 9; 38 15; 38 26; 39 7; 39 8; 39 11]
global d_x = [7.0, 5.0, 4.0, 3.0, 8.0, 8.0, 8.0, 9.0, 5.0, 8.0, 5.0, 7.0, 7.0, 1.0, 7.0, 5.0, 1.0, 9.0, 5.0, 6.0, 9.0, 1.0, 10.0, 7.0, 6.0, 5.0, 5.0, 7.0, 9.0, 4.0, 5.0, 8.0, 3.0, 7.0, 5.0, 10.0, 2.0, 7.0, 1.0, 1.0, 3.0, 1.0, 4.0, 6.0, 7.0, 2.0, 2.0, 3.0, 7.0, 9.0, 8.0, 8.0, 3.0, 4.0, 8.0, 8.0, 7.0, 9.0, 1.0, 6.0, 8.0, 7.0, 3.0, 2.0, 2.0, 8.0, 3.0, 7.0, 7.0, 9.0, 7.0, 2.0, 6.0, 10.0, 8.0, 4.0, 8.0, 10.0, 6.0, 6.0, 2.0, 7.0, 6.0, 3.0, 3.0, 2.0, 9.0, 9.0, 9.0, 5.0, 5.0, 9.0, 4.0, 5.0, 5.0, 7.0, 1.0, 3.0, 6.0, 1.0, 10.0, 9.0, 6.0, 7.0, 3.0, 1.0, 3.0, 10.0, 2.0, 4.0, 4.0, 5.0, 7.0, 3.0, 1.0, 5.0, 4.0, 10.0, 4.0, 3.0, 2.0, 2.0, 7.0, 6.0, 8.0, 7.0, 9.0, 2.0, 5.0, 8.0, 4.0, 10.0, 2.0, 2.0, 6.0, 2.0, 2.0, 2.0, 8.0, 1.0, 3.0, 10.0, 6.0, 5.0, 2.0, 6.0, 9.0, 1.0, 4.0, 3.0, 10.0, 7.0, 7.0, 8.0, 3.0, 3.0, 9.0, 4.0, 2.0, 2.0, 4.0, 6.0, 6.0, 7.0, 8.0, 10.0, 7.0, 5.0, 7.0]
global b_x = 5
global d_y = [8.0, 7.0, 1.0, 8.0, 4.0, 9.0, 5.0, 3.0, 2.0, 10.0, 9.0, 7.0, 3.0, 10.0, 3.0, 9.0, 4.0, 5.0, 1.0, 9.0, 8.0, 10.0, 2.0, 3.0, 4.0, 9.0, 6.0, 8.0, 5.0, 4.0, 2.0, 6.0, 5.0, 7.0, 2.0, 10.0, 6.0, 9.0, 9.0, 5.0, 9.0, 9.0, 6.0, 4.0, 4.0, 6.0, 5.0, 10.0, 7.0, 10.0, 3.0, 9.0, 1.0, 9.0, 5.0, 4.0, 2.0, 5.0, 1.0, 3.0, 8.0, 9.0, 6.0, 5.0, 6.0, 10.0, 1.0, 6.0, 5.0, 4.0, 9.0, 6.0, 7.0, 3.0, 7.0, 10.0, 4.0, 10.0, 1.0, 1.0, 9.0, 7.0, 3.0, 10.0, 2.0, 7.0, 10.0, 2.0, 7.0, 10.0, 8.0, 2.0, 1.0, 2.0, 8.0, 7.0, 3.0, 4.0, 2.0, 7.0, 3.0, 4.0, 9.0, 10.0, 6.0, 1.0, 8.0, 10.0, 4.0, 4.0, 7.0, 2.0, 4.0, 9.0, 1.0, 8.0, 3.0, 3.0, 10.0, 8.0, 5.0, 10.0, 10.0, 9.0, 2.0, 8.0, 7.0, 8.0, 1.0, 1.0, 3.0, 7.0, 10.0, 1.0, 8.0, 5.0, 9.0, 1.0, 3.0, 7.0, 9.0, 9.0, 6.0, 5.0, 4.0, 9.0, 2.0, 6.0, 8.0, 2.0, 6.0, 10.0, 6.0, 1.0, 6.0, 9.0, 4.0, 5.0, 2.0, 7.0, 2.0, 7.0, 5.0, 7.0, 4.0, 4.0, 10.0, 4.0, 5.0]
global b_y = 10
global p = [0.827, 0.472, 0.823, 0.626, 0.545, 0.246, 0.757, 0.428, 0.959, 0.302, 0.749, 0.664, 0.117, 0.527, 0.248, 0.35, 0.148, 0.726, 0.515, 0.411, 0.042, 0.181, 0.944, 0.456, 0.27, 0.949, 0.942, 0.223, 0.709, 0.177, 0.757, 0.094, 0.254, 0.141, 0.358, 0.307, 0.77, 0.522, 0.333, 0.017, 0.921, 0.34, 0.442, 0.095, 0.462, 0.036, 0.569, 0.852, 0.872, 0.096, 0.721, 0.656, 0.96, 0.814, 0.981, 0.331, 0.364, 0.483, 0.698, 0.5, 0.125, 0.708, 0.461, 0.055, 0.634, 0.03, 0.502, 0.806, 0.338, 0.673, 0.95, 0.638, 0.036, 0.932, 0.652, 0.572, 0.489, 0.357, 0.356, 0.002, 0.661, 0.465, 0.951, 0.308, 0.226, 0.143, 0.362, 0.038, 0.954, 0.661, 0.579, 0.975, 0.817, 0.063, 0.094, 0.919, 0.724, 0.892, 0.715, 0.786, 0.577, 0.621, 0.073, 0.128, 0.348, 0.934, 0.361, 0.003, 0.462, 0.913, 0.995, 0.482, 0.451, 0.689, 0.995, 0.783, 0.314, 0.455, 0.047, 0.767, 0.92, 0.614, 0.692, 0.411, 0.246, 0.302, 0.293, 0.73, 0.227, 0.823, 0.185, 0.799, 0.637, 0.344, 0.971, 0.42, 0.894, 0.721, 0.772, 0.871, 0.584, 0.536, 0.274, 0.898, 0.893, 0.246, 0.024, 0.387, 0.761, 0.543, 0.308, 0.144, 0.638, 0.889, 0.019, 0.403, 0.911, 0.879, 0.272, 0.21, 0.425, 0.482, 0.913, 0.588, 0.205, 0.537, 0.235, 0.736, 0.722]
global q = [0.967, 0.602, 0.894, 0.991, 0.865, 0.402, 0.91, 0.561, 0.978, 0.888, 0.905, 0.698, 0.914, 0.636, 0.362, 0.753, 0.968, 0.891, 0.934, 0.426, 0.635, 0.604, 0.978, 0.765, 0.716, 0.997, 0.958, 0.85, 0.963, 0.941, 0.978, 0.15, 0.369, 0.342, 0.484, 0.373, 0.983, 0.591, 0.969, 0.544, 0.964, 0.825, 0.84, 0.896, 0.946, 0.286, 0.799, 0.922, 0.938, 0.123, 0.957, 0.821, 0.97, 0.981, 0.99, 0.692, 0.627, 0.734, 0.983, 0.751, 0.287, 0.899, 0.48, 0.81, 0.907, 0.969, 0.657, 0.949, 0.711, 0.709, 0.963, 0.724, 0.728, 0.989, 0.668, 0.975, 0.556, 0.792, 0.923, 0.846, 0.718, 0.472, 0.961, 0.841, 0.688, 0.569, 0.484, 0.123, 0.993, 0.797, 0.584, 0.982, 0.934, 0.798, 0.846, 0.951, 0.733, 0.915, 0.931, 0.84, 0.588, 0.979, 0.905, 0.308, 0.604, 0.983, 0.46, 0.718, 0.61, 0.941, 0.995, 0.74, 0.639, 0.83, 0.995, 0.798, 0.468, 0.697, 0.517, 0.813, 0.947, 0.65, 0.987, 0.462, 0.354, 0.427, 0.902, 0.923, 0.907, 0.895, 0.892, 0.805, 0.755, 0.732, 0.972, 0.803, 0.912, 0.8, 0.929, 0.944, 0.625, 0.854, 0.56, 0.912, 0.963, 0.833, 0.083, 0.721, 0.999, 0.804, 0.341, 0.342, 0.725, 0.939, 0.301, 0.606, 0.925, 0.88, 0.308, 0.289, 0.892, 0.924, 0.967, 0.749, 0.641, 0.708, 0.537, 0.903, 0.739]
global origin = 1
global destination = 40