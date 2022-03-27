global arcs = [1 2; 1 29; 1 32; 1 35; 1 44; 1 46; 1 50; 2 10; 2 14; 2 16; 2 29; 2 36; 3 13; 3 21; 3 30; 3 38; 3 39; 3 40; 3 50; 4 3; 4 5; 4 6; 4 13; 4 23; 4 26; 4 34; 4 35; 4 36; 5 17; 5 23; 5 24; 5 35; 5 44; 6 14; 6 19; 6 25; 6 26; 6 40; 6 50; 7 19; 7 22; 7 23; 7 24; 7 25; 7 28; 7 30; 7 31; 7 45; 8 10; 8 43; 9 11; 9 18; 9 20; 9 30; 9 46; 10 11; 10 14; 10 16; 10 17; 10 21; 10 25; 10 39; 10 46; 10 50; 11 34; 11 39; 11 45; 11 50; 12 16; 13 12; 13 31; 13 33; 13 47; 14 15; 14 29; 14 37; 14 40; 15 7; 15 24; 15 25; 15 35; 15 39; 15 45; 16 5; 16 8; 16 20; 16 36; 16 44; 17 8; 17 27; 17 38; 17 45; 18 5; 18 14; 18 48; 19 11; 19 15; 19 27; 20 9; 20 12; 20 13; 20 46; 21 5; 21 7; 21 12; 21 14; 21 27; 21 33; 21 37; 21 38; 22 6; 22 10; 22 12; 22 25; 22 28; 22 41; 23 2; 23 4; 23 6; 23 16; 23 49; 24 19; 24 32; 24 50; 25 2; 25 13; 25 29; 25 48; 26 12; 26 22; 26 45; 27 3; 27 35; 27 39; 27 44; 28 29; 29 5; 29 20; 29 28; 29 40; 29 44; 30 2; 30 3; 30 13; 30 17; 30 21; 30 33; 30 39; 30 40; 31 4; 31 39; 31 45; 32 7; 32 31; 32 44; 33 4; 33 30; 33 37; 33 42; 34 25; 34 43; 34 45; 34 47; 35 42; 35 45; 36 30; 36 35; 36 41; 36 42; 36 49; 37 9; 37 19; 37 22; 37 28; 37 31; 37 39; 37 48; 38 9; 38 23; 38 34; 38 49; 38 50; 39 26; 39 35; 39 43; 40 5; 40 13; 40 20; 40 43; 41 5; 41 7; 41 12; 41 20; 41 39; 41 43; 41 50; 42 23; 42 30; 42 47; 43 4; 43 14; 43 24; 43 42; 43 47; 43 48; 44 2; 44 17; 44 20; 45 3; 45 11; 45 15; 45 23; 45 34; 46 17; 46 31; 46 39; 47 9; 47 13; 47 30; 47 35; 47 38; 47 48; 48 17; 48 41; 48 44; 48 45; 49 6; 49 9; 49 14; 49 16; 49 21; 49 28; 49 36]
global d_x = [8.0, 9.0, 2.0, 9.0, 2.0, 10.0, 1.0, 8.0, 7.0, 5.0, 3.0, 3.0, 5.0, 5.0, 7.0, 9.0, 10.0, 6.0, 6.0, 6.0, 8.0, 4.0, 10.0, 5.0, 9.0, 5.0, 4.0, 10.0, 5.0, 8.0, 2.0, 10.0, 4.0, 8.0, 4.0, 9.0, 4.0, 8.0, 1.0, 2.0, 2.0, 8.0, 5.0, 3.0, 7.0, 1.0, 7.0, 5.0, 2.0, 8.0, 6.0, 4.0, 9.0, 9.0, 3.0, 8.0, 5.0, 6.0, 4.0, 1.0, 3.0, 10.0, 2.0, 3.0, 7.0, 1.0, 4.0, 6.0, 8.0, 1.0, 9.0, 8.0, 9.0, 8.0, 2.0, 6.0, 1.0, 1.0, 9.0, 1.0, 6.0, 9.0, 7.0, 10.0, 9.0, 10.0, 3.0, 3.0, 4.0, 2.0, 1.0, 1.0, 10.0, 1.0, 3.0, 6.0, 10.0, 4.0, 8.0, 4.0, 5.0, 8.0, 10.0, 2.0, 5.0, 8.0, 2.0, 4.0, 10.0, 3.0, 7.0, 7.0, 9.0, 7.0, 6.0, 9.0, 3.0, 10.0, 2.0, 9.0, 3.0, 7.0, 9.0, 9.0, 7.0, 6.0, 4.0, 2.0, 2.0, 6.0, 1.0, 4.0, 3.0, 1.0, 5.0, 9.0, 5.0, 10.0, 1.0, 9.0, 9.0, 9.0, 1.0, 8.0, 8.0, 6.0, 3.0, 10.0, 1.0, 9.0, 10.0, 4.0, 8.0, 7.0, 2.0, 1.0, 10.0, 9.0, 6.0, 1.0, 5.0, 9.0, 4.0, 8.0, 5.0, 1.0, 1.0, 1.0, 10.0, 4.0, 6.0, 9.0, 7.0, 10.0, 9.0, 6.0, 4.0, 6.0, 5.0, 9.0, 7.0, 3.0, 7.0, 7.0, 10.0, 8.0, 7.0, 1.0, 5.0, 3.0, 9.0, 10.0, 6.0, 10.0, 1.0, 10.0, 8.0, 7.0, 9.0, 3.0, 2.0, 7.0, 9.0, 2.0, 10.0, 7.0, 5.0, 4.0, 1.0, 2.0, 10.0, 2.0, 6.0, 8.0, 5.0, 5.0, 3.0, 1.0, 10.0, 4.0, 7.0, 10.0, 10.0, 10.0, 1.0, 1.0, 3.0, 1.0, 6.0, 3.0, 8.0, 10.0, 6.0]
global b_x = 5
global d_y = [3.0, 10.0, 10.0, 8.0, 8.0, 5.0, 9.0, 6.0, 4.0, 5.0, 3.0, 1.0, 10.0, 5.0, 3.0, 4.0, 7.0, 10.0, 1.0, 5.0, 2.0, 9.0, 5.0, 9.0, 5.0, 2.0, 1.0, 2.0, 4.0, 10.0, 1.0, 3.0, 1.0, 8.0, 10.0, 4.0, 4.0, 3.0, 9.0, 3.0, 5.0, 8.0, 1.0, 7.0, 6.0, 8.0, 5.0, 5.0, 1.0, 7.0, 6.0, 7.0, 7.0, 10.0, 2.0, 10.0, 9.0, 6.0, 5.0, 2.0, 7.0, 8.0, 2.0, 3.0, 9.0, 2.0, 3.0, 10.0, 1.0, 5.0, 8.0, 6.0, 9.0, 8.0, 4.0, 10.0, 1.0, 2.0, 8.0, 5.0, 4.0, 3.0, 7.0, 10.0, 6.0, 10.0, 4.0, 3.0, 1.0, 3.0, 3.0, 4.0, 6.0, 8.0, 3.0, 9.0, 6.0, 7.0, 3.0, 2.0, 6.0, 3.0, 3.0, 10.0, 6.0, 3.0, 6.0, 7.0, 6.0, 7.0, 7.0, 9.0, 2.0, 9.0, 3.0, 9.0, 7.0, 4.0, 6.0, 5.0, 5.0, 8.0, 1.0, 10.0, 7.0, 7.0, 1.0, 10.0, 1.0, 9.0, 10.0, 3.0, 8.0, 7.0, 7.0, 6.0, 1.0, 3.0, 9.0, 7.0, 5.0, 4.0, 7.0, 8.0, 6.0, 9.0, 8.0, 3.0, 3.0, 2.0, 1.0, 3.0, 9.0, 4.0, 1.0, 8.0, 10.0, 4.0, 1.0, 1.0, 9.0, 4.0, 10.0, 7.0, 7.0, 5.0, 6.0, 8.0, 8.0, 6.0, 4.0, 4.0, 2.0, 7.0, 8.0, 2.0, 10.0, 9.0, 6.0, 10.0, 5.0, 2.0, 8.0, 1.0, 9.0, 5.0, 8.0, 4.0, 8.0, 3.0, 8.0, 5.0, 9.0, 9.0, 6.0, 1.0, 7.0, 9.0, 5.0, 7.0, 7.0, 3.0, 7.0, 3.0, 6.0, 3.0, 8.0, 2.0, 7.0, 10.0, 8.0, 5.0, 8.0, 4.0, 4.0, 10.0, 3.0, 8.0, 8.0, 2.0, 10.0, 2.0, 4.0, 2.0, 7.0, 7.0, 4.0, 2.0, 4.0, 6.0, 3.0, 8.0, 1.0]
global b_y = 10
global p = [0.581, 0.568, 0.061, 0.068, 0.686, 0.155, 0.371, 0.384, 0.571, 0.223, 0.353, 0.777, 0.268, 0.616, 0.826, 0.134, 0.701, 0.551, 0.918, 0.474, 0.196, 0.395, 0.999, 0.405, 0.931, 0.145, 0.48, 0.726, 0.347, 0.546, 0.894, 0.23, 0.454, 0.89, 0.131, 0.444, 0.187, 0.744, 0.263, 0.421, 0.103, 0.787, 0.673, 0.608, 0.811, 0.832, 0.591, 0.768, 0.572, 0.084, 0.512, 0.398, 0.055, 0.031, 0.412, 0.046, 0.728, 0.979, 0.447, 0.745, 0.307, 0.235, 0.323, 0.827, 0.187, 0.334, 0.937, 0.281, 0.005, 0.744, 0.218, 0.942, 0.167, 0.282, 0.383, 0.835, 0.605, 0.503, 0.169, 0.851, 0.061, 0.007, 0.732, 0.855, 0.618, 0.073, 0.799, 0.026, 0.322, 0.605, 0.019, 0.573, 0.001, 0.244, 0.144, 0.754, 0.924, 0.072, 0.456, 0.439, 0.749, 0.728, 0.18, 0.127, 0.629, 0.048, 0.486, 0.013, 0.784, 0.846, 0.398, 0.306, 0.007, 0.787, 0.948, 0.389, 0.552, 0.633, 0.941, 0.759, 0.794, 0.358, 0.124, 0.07, 0.848, 0.418, 0.916, 0.074, 0.533, 0.635, 0.958, 0.837, 0.631, 0.261, 0.629, 0.931, 0.76, 0.009, 0.157, 0.764, 0.137, 0.293, 0.516, 0.685, 0.045, 0.642, 0.852, 0.398, 0.036, 0.615, 0.767, 0.59, 0.388, 0.722, 0.502, 0.91, 0.595, 0.107, 0.014, 0.087, 0.737, 0.111, 0.709, 0.984, 0.326, 0.942, 0.462, 0.38, 0.254, 0.386, 0.313, 0.215, 0.498, 0.927, 0.467, 0.221, 0.378, 0.72, 0.615, 0.806, 0.962, 0.073, 0.952, 0.117, 0.167, 0.54, 0.353, 0.987, 0.714, 0.168, 0.397, 0.285, 0.264, 0.072, 0.783, 0.485, 0.229, 0.101, 0.512, 0.003, 0.557, 0.729, 0.699, 0.322, 0.857, 0.581, 0.889, 0.045, 0.22, 0.096, 0.298, 0.143, 0.075, 0.382, 0.182, 0.639, 0.332, 0.464, 0.543, 0.654, 0.221, 0.704, 0.865, 0.409, 0.102, 0.379, 0.27, 0.24, 0.14, 0.624, 0.84, 0.027, 0.9]
global q = [0.99, 0.595, 0.425, 0.976, 0.703, 0.68, 0.641, 0.613, 0.848, 0.716, 0.977, 0.874, 0.806, 0.771, 0.934, 0.624, 0.736, 0.741, 0.991, 0.484, 0.452, 0.814, 0.999, 0.785, 0.993, 0.873, 0.625, 0.819, 0.691, 0.914, 0.977, 0.503, 0.893, 0.941, 0.244, 0.679, 0.955, 0.887, 0.798, 0.997, 0.678, 0.951, 0.747, 0.823, 0.9, 0.848, 0.809, 0.949, 0.947, 0.929, 0.751, 0.556, 0.858, 0.774, 0.762, 0.701, 0.97, 0.994, 0.82, 0.8, 0.458, 0.883, 0.52, 0.947, 0.205, 0.494, 0.957, 0.481, 0.255, 0.849, 0.745, 0.951, 0.374, 0.891, 0.84, 0.947, 0.88, 0.552, 0.579, 0.922, 0.8, 0.658, 0.99, 0.873, 0.816, 0.147, 0.934, 0.387, 0.362, 0.969, 0.117, 0.841, 0.344, 0.684, 0.48, 0.783, 0.947, 0.265, 0.477, 0.683, 0.848, 0.744, 0.427, 0.952, 0.807, 0.349, 0.675, 0.696, 0.793, 0.991, 0.842, 0.702, 0.286, 0.992, 0.95, 0.952, 0.775, 0.837, 0.963, 0.815, 0.916, 0.528, 0.81, 0.401, 0.883, 0.965, 0.967, 0.433, 0.82, 0.999, 0.995, 0.864, 0.763, 0.689, 0.664, 0.98, 0.947, 0.345, 0.461, 0.862, 0.75, 0.714, 0.864, 0.705, 0.721, 0.935, 0.893, 0.804, 0.953, 0.62, 0.977, 0.869, 0.556, 0.898, 0.811, 0.989, 0.651, 0.362, 0.284, 0.851, 0.786, 0.961, 0.774, 0.998, 0.606, 0.957, 0.973, 0.711, 0.569, 0.388, 0.544, 0.993, 0.947, 0.956, 0.992, 0.711, 0.667, 0.795, 0.84, 0.885, 0.97, 0.706, 0.955, 0.563, 0.577, 0.553, 0.881, 0.987, 0.793, 0.81, 0.397, 0.285, 0.839, 0.743, 0.858, 0.489, 0.5, 0.78, 0.73, 0.767, 0.876, 0.837, 0.944, 0.548, 0.932, 0.906, 0.969, 0.482, 0.928, 0.611, 0.419, 0.719, 0.942, 0.863, 0.194, 0.852, 0.588, 0.926, 0.92, 0.81, 0.997, 0.893, 0.907, 0.806, 0.69, 0.906, 0.389, 0.522, 0.494, 0.648, 0.978, 0.269, 0.914]
global origin = 1
global destination = 50