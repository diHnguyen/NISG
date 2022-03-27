global arcs = [1 5; 1 15; 1 19; 1 22; 1 39; 1 40; 1 45; 2 12; 3 6; 3 11; 3 19; 3 29; 3 34; 4 15; 4 32; 4 34; 4 36; 4 45; 5 6; 5 14; 5 25; 5 37; 5 49; 6 8; 6 10; 6 24; 6 27; 6 49; 7 5; 7 9; 7 11; 7 13; 7 26; 7 27; 7 38; 8 9; 8 20; 8 42; 8 43; 9 21; 9 24; 9 26; 9 29; 9 34; 9 39; 9 48; 10 6; 10 16; 10 21; 10 36; 10 43; 11 3; 11 16; 11 23; 11 30; 11 31; 11 38; 11 43; 12 10; 12 13; 12 16; 12 21; 12 46; 12 47; 12 49; 12 50; 13 11; 13 12; 13 28; 13 38; 13 40; 14 11; 14 12; 14 26; 14 27; 14 28; 14 33; 15 10; 15 43; 16 5; 16 10; 16 13; 16 18; 16 32; 16 35; 16 42; 17 2; 17 4; 17 18; 17 33; 17 34; 17 36; 17 37; 18 2; 18 5; 18 22; 18 26; 18 27; 18 36; 18 43; 18 44; 18 45; 18 46; 19 2; 19 28; 19 29; 19 30; 19 44; 19 46; 20 17; 20 28; 20 31; 20 35; 20 50; 21 14; 21 17; 21 38; 21 47; 22 6; 22 10; 22 34; 23 4; 23 16; 23 19; 23 28; 23 36; 23 41; 24 4; 25 24; 25 32; 25 44; 25 49; 26 22; 26 35; 27 10; 27 22; 27 33; 28 4; 28 13; 28 21; 28 39; 28 42; 28 45; 29 16; 29 17; 29 22; 30 18; 30 25; 30 33; 30 34; 30 39; 31 15; 31 19; 32 8; 32 24; 32 31; 32 36; 32 37; 32 38; 32 39; 32 44; 32 45; 33 11; 33 20; 33 22; 33 31; 33 41; 34 5; 34 20; 34 24; 34 31; 34 48; 35 10; 35 20; 35 23; 35 36; 35 44; 35 48; 36 9; 36 22; 36 30; 36 44; 36 50; 37 9; 37 33; 37 43; 38 24; 38 26; 38 30; 39 15; 39 32; 39 33; 40 7; 40 8; 40 25; 40 29; 40 34; 40 42; 40 48; 40 49; 41 4; 41 15; 41 24; 42 6; 42 9; 43 16; 43 21; 43 25; 43 39; 43 46; 44 10; 44 22; 44 32; 44 35; 45 2; 45 12; 45 17; 45 29; 45 31; 45 33; 45 48; 46 7; 46 20; 46 21; 46 27; 46 34; 46 39; 47 8; 47 12; 47 22; 47 25; 47 38; 48 5; 48 15; 48 18; 48 29; 48 30; 48 32; 48 36; 48 40; 48 42; 48 44; 49 3; 49 10; 49 14; 49 18; 49 29; 49 34; 49 41; 49 50]
global d_x = [3.0, 10.0, 9.0, 9.0, 8.0, 4.0, 10.0, 5.0, 9.0, 8.0, 3.0, 7.0, 2.0, 9.0, 6.0, 6.0, 2.0, 2.0, 3.0, 8.0, 4.0, 3.0, 10.0, 2.0, 8.0, 2.0, 10.0, 4.0, 2.0, 9.0, 8.0, 4.0, 9.0, 3.0, 3.0, 4.0, 10.0, 5.0, 9.0, 4.0, 2.0, 5.0, 2.0, 2.0, 3.0, 7.0, 2.0, 10.0, 7.0, 8.0, 5.0, 7.0, 10.0, 10.0, 8.0, 4.0, 7.0, 7.0, 1.0, 3.0, 2.0, 9.0, 6.0, 5.0, 10.0, 10.0, 4.0, 6.0, 10.0, 9.0, 8.0, 10.0, 10.0, 2.0, 2.0, 6.0, 5.0, 9.0, 6.0, 2.0, 10.0, 6.0, 9.0, 9.0, 4.0, 9.0, 2.0, 6.0, 1.0, 2.0, 8.0, 10.0, 4.0, 6.0, 5.0, 6.0, 4.0, 5.0, 9.0, 8.0, 10.0, 4.0, 5.0, 10.0, 6.0, 7.0, 5.0, 9.0, 6.0, 8.0, 4.0, 2.0, 8.0, 5.0, 7.0, 10.0, 9.0, 8.0, 2.0, 7.0, 4.0, 5.0, 10.0, 1.0, 4.0, 3.0, 1.0, 3.0, 6.0, 7.0, 5.0, 4.0, 5.0, 4.0, 7.0, 6.0, 2.0, 7.0, 5.0, 8.0, 8.0, 1.0, 1.0, 9.0, 2.0, 9.0, 8.0, 3.0, 10.0, 3.0, 7.0, 7.0, 8.0, 9.0, 5.0, 10.0, 6.0, 5.0, 1.0, 5.0, 7.0, 4.0, 3.0, 6.0, 6.0, 1.0, 2.0, 7.0, 6.0, 8.0, 7.0, 3.0, 6.0, 9.0, 5.0, 1.0, 3.0, 2.0, 5.0, 7.0, 3.0, 8.0, 7.0, 6.0, 4.0, 5.0, 6.0, 5.0, 2.0, 10.0, 5.0, 6.0, 8.0, 7.0, 8.0, 9.0, 6.0, 4.0, 10.0, 9.0, 9.0, 10.0, 1.0, 8.0, 3.0, 4.0, 7.0, 8.0, 5.0, 9.0, 6.0, 1.0, 3.0, 1.0, 9.0, 6.0, 6.0, 6.0, 5.0, 2.0, 2.0, 3.0, 6.0, 4.0, 5.0, 5.0, 4.0, 6.0, 1.0, 6.0, 10.0, 4.0, 10.0, 5.0, 7.0, 8.0, 5.0, 4.0, 3.0, 1.0, 6.0, 4.0, 7.0, 6.0, 6.0, 9.0, 9.0, 2.0, 8.0, 3.0]
global b_x = 5
global d_y = [10.0, 4.0, 10.0, 9.0, 1.0, 8.0, 6.0, 9.0, 8.0, 9.0, 1.0, 7.0, 9.0, 6.0, 10.0, 10.0, 1.0, 8.0, 2.0, 8.0, 2.0, 8.0, 2.0, 3.0, 2.0, 2.0, 1.0, 9.0, 7.0, 5.0, 8.0, 2.0, 8.0, 5.0, 5.0, 9.0, 10.0, 10.0, 5.0, 8.0, 5.0, 6.0, 8.0, 2.0, 5.0, 4.0, 5.0, 8.0, 2.0, 3.0, 5.0, 5.0, 2.0, 4.0, 7.0, 7.0, 2.0, 7.0, 7.0, 4.0, 6.0, 4.0, 9.0, 7.0, 5.0, 6.0, 6.0, 1.0, 8.0, 6.0, 10.0, 3.0, 5.0, 10.0, 3.0, 8.0, 9.0, 8.0, 10.0, 3.0, 10.0, 4.0, 2.0, 6.0, 5.0, 1.0, 5.0, 3.0, 4.0, 1.0, 5.0, 7.0, 8.0, 7.0, 10.0, 10.0, 4.0, 2.0, 1.0, 9.0, 8.0, 9.0, 6.0, 1.0, 3.0, 7.0, 3.0, 6.0, 1.0, 8.0, 8.0, 7.0, 10.0, 9.0, 7.0, 8.0, 1.0, 3.0, 5.0, 8.0, 9.0, 10.0, 4.0, 10.0, 10.0, 6.0, 2.0, 4.0, 6.0, 10.0, 3.0, 3.0, 2.0, 6.0, 10.0, 4.0, 7.0, 8.0, 4.0, 2.0, 4.0, 2.0, 4.0, 10.0, 6.0, 5.0, 7.0, 6.0, 10.0, 2.0, 1.0, 8.0, 8.0, 8.0, 5.0, 9.0, 10.0, 8.0, 3.0, 7.0, 6.0, 2.0, 8.0, 8.0, 4.0, 6.0, 8.0, 1.0, 4.0, 9.0, 5.0, 4.0, 5.0, 2.0, 2.0, 3.0, 1.0, 1.0, 6.0, 3.0, 7.0, 9.0, 3.0, 4.0, 1.0, 6.0, 10.0, 9.0, 1.0, 6.0, 2.0, 5.0, 6.0, 3.0, 2.0, 7.0, 9.0, 7.0, 3.0, 7.0, 2.0, 4.0, 4.0, 2.0, 1.0, 7.0, 9.0, 10.0, 1.0, 6.0, 6.0, 8.0, 6.0, 3.0, 4.0, 5.0, 1.0, 3.0, 10.0, 5.0, 10.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 10.0, 5.0, 8.0, 2.0, 4.0, 9.0, 5.0, 4.0, 7.0, 3.0, 1.0, 10.0, 4.0, 1.0, 6.0, 7.0, 9.0, 8.0, 9.0, 5.0, 6.0, 2.0, 6.0]
global b_y = 10
global p = [0.458, 0.04, 0.105, 0.884, 0.35, 0.448, 0.543, 0.897, 0.035, 0.002, 0.613, 0.021, 0.556, 0.562, 0.131, 0.424, 0.541, 0.599, 0.961, 0.794, 0.278, 0.436, 0.56, 0.107, 0.161, 0.987, 0.98, 0.782, 0.051, 0.618, 0.37, 0.883, 0.371, 0.982, 0.557, 0.305, 0.465, 0.829, 0.843, 0.29, 0.693, 0.649, 0.177, 0.437, 0.461, 0.548, 0.132, 0.53, 0.641, 0.56, 0.787, 0.43, 0.738, 0.783, 0.335, 0.594, 0.514, 0.24, 0.165, 0.028, 0.65, 0.828, 0.018, 0.209, 0.013, 0.147, 0.486, 0.805, 0.767, 0.733, 0.962, 0.297, 0.051, 0.326, 0.11, 0.696, 0.566, 0.918, 0.985, 0.259, 0.226, 0.983, 0.9, 0.951, 0.658, 0.965, 0.451, 0.61, 0.859, 0.784, 0.72, 0.972, 0.388, 0.225, 0.011, 0.632, 0.536, 0.749, 0.755, 0.172, 0.675, 0.729, 0.944, 0.785, 0.133, 0.706, 0.573, 0.266, 0.357, 0.47, 0.683, 0.871, 0.245, 0.261, 0.38, 0.408, 0.536, 0.813, 0.486, 0.846, 0.322, 0.095, 0.667, 0.469, 0.104, 0.76, 0.187, 0.871, 0.283, 0.619, 0.499, 0.298, 0.015, 0.799, 0.067, 0.509, 0.959, 0.412, 0.939, 0.734, 0.742, 0.777, 0.231, 0.152, 0.453, 0.41, 0.522, 0.497, 0.243, 0.446, 0.556, 0.587, 0.01, 0.183, 0.005, 0.013, 0.599, 0.217, 0.422, 0.373, 0.166, 0.609, 0.863, 0.586, 0.583, 0.729, 0.833, 0.598, 0.274, 0.755, 0.671, 0.859, 0.152, 0.026, 0.55, 0.605, 0.039, 0.926, 0.942, 0.44, 0.813, 0.412, 0.538, 0.179, 0.23, 0.564, 0.876, 0.243, 0.473, 0.655, 0.22, 0.027, 0.042, 0.411, 0.334, 0.912, 0.831, 0.374, 0.242, 0.836, 0.895, 0.399, 0.105, 0.245, 0.354, 0.218, 0.672, 0.128, 0.645, 0.598, 0.357, 0.185, 0.937, 0.709, 0.371, 0.471, 0.745, 0.356, 0.508, 0.412, 0.81, 0.809, 0.519, 0.09, 0.134, 0.181, 0.926, 0.542, 0.477, 0.914, 0.127, 0.182, 0.796, 0.504, 0.762, 0.796, 0.23, 0.874, 0.801, 0.415, 0.84, 0.524, 0.256, 0.892, 0.981, 0.075, 0.282, 0.163, 0.836, 0.767]
global q = [0.767, 0.146, 0.116, 0.963, 0.576, 0.563, 0.862, 0.952, 0.473, 0.318, 0.964, 0.446, 0.613, 0.597, 0.938, 0.669, 0.871, 0.955, 0.982, 0.926, 0.747, 0.942, 0.786, 0.931, 0.177, 0.996, 0.997, 0.849, 0.117, 0.711, 0.497, 0.981, 0.8, 0.999, 0.882, 0.952, 0.54, 0.849, 0.944, 0.29, 0.844, 0.974, 0.587, 0.97, 0.693, 0.941, 0.843, 0.876, 0.739, 0.766, 0.898, 0.75, 0.981, 0.875, 0.661, 0.827, 0.549, 0.565, 0.425, 0.082, 0.762, 0.919, 0.688, 0.524, 0.604, 0.899, 0.517, 0.808, 0.908, 0.773, 0.969, 0.384, 0.847, 0.457, 0.626, 0.998, 0.668, 0.974, 0.995, 0.684, 0.748, 0.995, 0.917, 0.987, 0.879, 0.967, 0.878, 0.687, 0.904, 0.819, 0.994, 0.984, 0.394, 0.341, 0.984, 0.801, 0.885, 0.883, 0.836, 0.575, 0.731, 0.814, 0.957, 0.914, 0.381, 0.781, 0.936, 0.364, 0.553, 0.734, 0.972, 0.971, 0.711, 0.909, 0.502, 0.525, 0.871, 0.848, 0.96, 0.931, 0.774, 0.242, 0.822, 0.921, 0.159, 0.913, 0.262, 0.908, 0.381, 0.659, 0.845, 0.434, 0.83, 0.805, 0.561, 0.943, 0.964, 0.658, 0.942, 0.911, 0.827, 0.796, 0.49, 0.437, 0.594, 0.487, 0.795, 0.503, 0.302, 0.649, 0.726, 0.868, 0.133, 0.52, 0.138, 0.749, 0.802, 0.77, 0.981, 0.41, 0.784, 0.643, 0.943, 0.721, 0.809, 0.898, 0.842, 0.877, 0.61, 0.828, 0.836, 0.906, 0.295, 0.521, 0.927, 0.658, 0.52, 0.989, 0.944, 0.605, 0.832, 0.694, 0.61, 0.573, 0.796, 0.788, 0.919, 0.729, 0.939, 0.688, 0.801, 0.338, 0.326, 0.818, 0.884, 0.968, 0.878, 0.468, 0.647, 0.924, 0.9, 0.723, 0.126, 0.784, 0.878, 0.786, 0.752, 0.663, 0.833, 0.849, 0.907, 0.871, 0.948, 0.848, 0.959, 0.651, 0.837, 0.739, 0.807, 0.468, 0.872, 0.939, 0.648, 0.352, 0.964, 0.63, 0.946, 0.752, 0.518, 0.92, 0.601, 0.474, 0.89, 0.804, 0.886, 0.97, 0.503, 0.968, 0.831, 0.807, 0.931, 0.53, 0.274, 0.969, 0.998, 0.147, 0.63, 0.319, 0.845, 0.921]
global origin = 1
global destination = 50