global arcs = [1 16; 1 32; 2 12; 2 16; 2 29; 3 6; 3 17; 3 43; 4 12; 4 13; 4 17; 4 49; 5 3; 5 11; 5 28; 5 29; 5 31; 5 34; 5 45; 6 10; 6 19; 6 28; 6 44; 7 2; 7 5; 7 15; 7 20; 7 39; 8 4; 8 11; 8 16; 8 20; 8 24; 8 30; 8 33; 8 35; 8 40; 9 13; 9 27; 9 28; 9 45; 9 49; 10 6; 10 42; 10 45; 10 46; 10 49; 11 2; 11 3; 11 4; 11 10; 11 34; 11 40; 11 48; 12 11; 12 18; 12 21; 12 29; 13 8; 13 24; 13 37; 14 4; 14 11; 14 12; 14 25; 14 37; 15 9; 15 14; 15 24; 15 25; 15 32; 15 35; 15 41; 16 19; 16 30; 16 40; 16 44; 16 49; 17 4; 17 15; 17 25; 17 30; 17 45; 18 27; 19 2; 19 8; 19 13; 19 15; 19 26; 19 36; 19 42; 20 17; 20 28; 20 32; 20 36; 20 42; 21 26; 21 29; 21 39; 21 42; 21 45; 21 47; 22 3; 22 18; 22 23; 22 31; 22 39; 23 7; 24 2; 24 4; 24 5; 24 14; 24 33; 25 8; 25 19; 25 41; 25 46; 26 2; 26 18; 26 19; 26 24; 26 29; 26 32; 27 6; 27 30; 27 34; 27 48; 28 8; 28 12; 28 13; 28 36; 28 42; 29 14; 29 19; 29 23; 29 32; 29 42; 30 17; 30 39; 30 44; 30 47; 30 50; 31 9; 31 21; 31 26; 31 42; 32 4; 32 20; 32 21; 32 22; 32 26; 32 30; 32 37; 32 38; 33 3; 33 21; 33 28; 33 37; 33 45; 34 3; 34 17; 34 19; 34 31; 35 3; 35 19; 35 34; 35 42; 36 4; 36 7; 36 18; 36 25; 36 30; 36 32; 36 34; 36 44; 36 47; 37 9; 38 5; 38 23; 38 37; 38 42; 38 46; 39 20; 39 21; 39 41; 39 46; 39 49; 40 12; 40 15; 40 18; 40 28; 40 39; 40 48; 41 2; 41 9; 41 10; 41 12; 41 13; 41 14; 41 15; 41 29; 41 39; 41 44; 42 7; 42 41; 43 46; 44 9; 44 17; 44 23; 44 43; 45 18; 45 36; 45 40; 45 41; 45 42; 45 50; 46 5; 46 28; 47 5; 47 15; 47 25; 47 31; 47 39; 48 19; 48 23; 48 24; 48 26; 48 34; 48 42; 48 45; 49 2; 49 11; 49 38; 49 42; 49 46]
global d_x = [5.0, 9.0, 1.0, 5.0, 4.0, 9.0, 7.0, 8.0, 8.0, 10.0, 7.0, 2.0, 5.0, 10.0, 8.0, 6.0, 6.0, 6.0, 2.0, 10.0, 6.0, 5.0, 2.0, 9.0, 8.0, 7.0, 10.0, 2.0, 2.0, 2.0, 7.0, 2.0, 3.0, 7.0, 1.0, 1.0, 3.0, 6.0, 9.0, 6.0, 8.0, 1.0, 7.0, 5.0, 4.0, 4.0, 1.0, 7.0, 1.0, 8.0, 1.0, 10.0, 2.0, 4.0, 10.0, 1.0, 3.0, 7.0, 7.0, 8.0, 9.0, 5.0, 3.0, 4.0, 5.0, 4.0, 9.0, 5.0, 5.0, 6.0, 10.0, 4.0, 6.0, 5.0, 6.0, 3.0, 10.0, 1.0, 6.0, 4.0, 2.0, 3.0, 6.0, 7.0, 3.0, 1.0, 5.0, 6.0, 10.0, 7.0, 10.0, 6.0, 2.0, 1.0, 2.0, 1.0, 3.0, 3.0, 8.0, 5.0, 5.0, 4.0, 3.0, 9.0, 6.0, 2.0, 10.0, 8.0, 9.0, 4.0, 6.0, 3.0, 8.0, 5.0, 10.0, 3.0, 3.0, 8.0, 4.0, 3.0, 2.0, 7.0, 3.0, 6.0, 8.0, 5.0, 7.0, 9.0, 2.0, 5.0, 2.0, 1.0, 4.0, 9.0, 8.0, 6.0, 8.0, 10.0, 3.0, 2.0, 8.0, 10.0, 3.0, 4.0, 5.0, 5.0, 3.0, 9.0, 2.0, 7.0, 7.0, 3.0, 1.0, 9.0, 7.0, 4.0, 6.0, 8.0, 5.0, 10.0, 5.0, 1.0, 1.0, 5.0, 1.0, 2.0, 1.0, 2.0, 5.0, 6.0, 9.0, 5.0, 9.0, 8.0, 9.0, 4.0, 1.0, 10.0, 4.0, 5.0, 7.0, 10.0, 1.0, 5.0, 5.0, 7.0, 1.0, 6.0, 9.0, 7.0, 5.0, 1.0, 3.0, 1.0, 2.0, 10.0, 8.0, 2.0, 8.0, 9.0, 4.0, 1.0, 4.0, 1.0, 6.0, 4.0, 1.0, 10.0, 6.0, 9.0, 10.0, 9.0, 4.0, 4.0, 9.0, 7.0, 8.0, 5.0, 10.0, 6.0, 6.0, 5.0, 1.0, 10.0, 4.0, 10.0, 3.0, 3.0, 7.0, 2.0, 5.0, 4.0, 6.0, 7.0, 8.0]
global b_x = 5
global d_y = [3.0, 7.0, 9.0, 8.0, 2.0, 9.0, 4.0, 5.0, 1.0, 3.0, 6.0, 7.0, 4.0, 10.0, 2.0, 2.0, 3.0, 10.0, 6.0, 6.0, 3.0, 3.0, 5.0, 7.0, 10.0, 6.0, 10.0, 2.0, 3.0, 6.0, 7.0, 2.0, 1.0, 8.0, 5.0, 4.0, 10.0, 7.0, 4.0, 4.0, 6.0, 4.0, 6.0, 3.0, 1.0, 1.0, 7.0, 2.0, 7.0, 7.0, 2.0, 2.0, 8.0, 8.0, 9.0, 8.0, 9.0, 2.0, 2.0, 10.0, 10.0, 2.0, 10.0, 8.0, 1.0, 6.0, 9.0, 4.0, 9.0, 3.0, 7.0, 3.0, 9.0, 10.0, 7.0, 10.0, 9.0, 5.0, 6.0, 3.0, 6.0, 2.0, 7.0, 2.0, 8.0, 10.0, 9.0, 10.0, 1.0, 10.0, 2.0, 6.0, 4.0, 1.0, 8.0, 4.0, 1.0, 5.0, 8.0, 4.0, 6.0, 8.0, 2.0, 1.0, 7.0, 6.0, 2.0, 1.0, 8.0, 6.0, 6.0, 1.0, 8.0, 10.0, 4.0, 8.0, 3.0, 1.0, 5.0, 1.0, 6.0, 9.0, 10.0, 3.0, 3.0, 6.0, 6.0, 7.0, 9.0, 6.0, 7.0, 1.0, 4.0, 9.0, 3.0, 1.0, 1.0, 10.0, 7.0, 5.0, 7.0, 1.0, 1.0, 2.0, 1.0, 8.0, 7.0, 3.0, 6.0, 6.0, 1.0, 7.0, 4.0, 10.0, 5.0, 6.0, 6.0, 2.0, 2.0, 3.0, 9.0, 5.0, 8.0, 5.0, 2.0, 7.0, 8.0, 10.0, 7.0, 2.0, 6.0, 6.0, 3.0, 8.0, 9.0, 9.0, 7.0, 6.0, 4.0, 9.0, 1.0, 3.0, 8.0, 9.0, 1.0, 6.0, 8.0, 7.0, 8.0, 4.0, 5.0, 6.0, 1.0, 8.0, 2.0, 6.0, 10.0, 6.0, 9.0, 6.0, 6.0, 3.0, 1.0, 1.0, 4.0, 1.0, 7.0, 3.0, 8.0, 1.0, 1.0, 8.0, 5.0, 6.0, 2.0, 6.0, 3.0, 9.0, 6.0, 3.0, 4.0, 7.0, 3.0, 9.0, 3.0, 10.0, 7.0, 6.0, 7.0, 3.0, 4.0, 1.0, 5.0, 7.0, 7.0]
global b_y = 10
global p = [0.883, 0.094, 0.757, 0.274, 0.733, 0.38, 0.412, 0.697, 0.042, 0.272, 0.709, 0.836, 0.023, 0.387, 0.098, 0.122, 0.598, 0.471, 0.482, 0.435, 0.585, 0.523, 0.206, 0.715, 0.051, 0.251, 0.706, 0.23, 0.038, 0.199, 0.04, 0.269, 0.848, 0.224, 0.768, 0.289, 0.608, 0.621, 0.012, 0.125, 0.61, 0.656, 0.447, 0.663, 0.55, 0.27, 0.599, 0.474, 0.129, 0.712, 0.163, 0.805, 0.471, 0.706, 0.067, 0.884, 0.889, 0.1, 0.5, 0.936, 0.859, 0.006, 0.205, 0.644, 0.215, 0.673, 0.075, 0.77, 0.778, 0.83, 0.93, 0.664, 0.551, 0.442, 0.291, 0.83, 0.169, 0.184, 0.965, 0.735, 0.921, 0.571, 0.041, 0.601, 0.225, 0.31, 0.081, 0.759, 0.099, 0.8, 0.821, 0.756, 0.563, 0.576, 0.242, 0.915, 0.871, 0.879, 0.82, 0.428, 0.386, 0.023, 0.886, 0.263, 0.946, 0.717, 0.234, 0.666, 0.902, 0.911, 0.3, 0.575, 0.655, 0.71, 0.464, 0.33, 0.85, 0.408, 0.294, 0.964, 0.407, 0.547, 0.983, 0.66, 0.077, 0.24, 0.686, 0.882, 0.211, 0.685, 0.125, 0.068, 0.069, 0.699, 0.889, 0.447, 0.109, 0.143, 0.883, 0.832, 0.399, 0.705, 0.609, 0.578, 0.254, 0.48, 0.019, 0.058, 0.975, 0.342, 0.962, 0.381, 0.693, 0.632, 0.027, 0.298, 0.255, 0.601, 0.501, 0.909, 0.635, 0.548, 0.394, 0.289, 0.03, 0.871, 0.594, 0.896, 0.775, 0.415, 0.932, 0.577, 0.03, 0.849, 0.455, 0.166, 0.731, 0.691, 0.36, 0.824, 0.887, 0.153, 0.828, 0.815, 0.616, 0.079, 0.124, 0.199, 0.965, 0.594, 0.269, 0.216, 0.084, 0.951, 0.799, 0.086, 0.777, 0.341, 0.318, 0.594, 0.359, 0.31, 0.296, 0.392, 0.476, 0.34, 0.662, 0.015, 0.283, 0.408, 0.292, 0.107, 0.075, 0.079, 0.544, 0.378, 0.493, 0.255, 0.186, 0.313, 0.308, 0.872, 0.357, 0.157, 0.823, 0.99, 0.148, 0.202, 0.116, 0.496, 0.183, 0.969, 0.089, 0.188, 0.392]
global q = [0.978, 0.975, 0.928, 0.871, 0.98, 0.571, 0.688, 0.99, 0.868, 0.935, 0.942, 0.938, 0.758, 0.65, 0.115, 0.768, 0.672, 0.918, 0.627, 0.668, 0.678, 0.653, 0.669, 0.728, 0.205, 0.587, 0.902, 0.731, 0.812, 0.924, 0.233, 0.645, 0.998, 0.68, 0.871, 0.688, 0.837, 0.756, 0.978, 0.764, 0.796, 0.969, 0.865, 0.773, 0.787, 0.78, 0.61, 0.934, 0.181, 0.956, 0.677, 0.977, 0.964, 0.998, 0.403, 0.903, 0.997, 0.55, 0.843, 0.996, 0.967, 0.402, 0.613, 0.661, 0.883, 0.934, 0.576, 0.898, 0.951, 0.977, 0.953, 0.745, 0.655, 0.592, 0.384, 0.916, 0.438, 0.722, 0.972, 0.793, 0.997, 0.737, 0.24, 0.633, 0.934, 0.776, 0.269, 0.965, 0.672, 0.975, 0.941, 0.971, 0.571, 0.576, 0.628, 0.973, 0.951, 0.922, 0.939, 0.691, 0.624, 0.215, 0.917, 0.697, 0.996, 0.722, 0.305, 0.983, 0.979, 0.957, 0.937, 0.714, 0.931, 0.764, 0.804, 0.845, 0.854, 0.816, 0.561, 0.981, 0.455, 0.743, 0.99, 0.708, 0.306, 0.653, 0.942, 0.934, 0.772, 0.699, 0.897, 0.438, 0.899, 0.778, 0.943, 0.967, 0.823, 0.746, 0.917, 0.883, 0.502, 0.908, 0.791, 0.733, 0.275, 0.835, 0.414, 0.661, 0.982, 0.372, 0.994, 0.648, 0.88, 0.769, 0.899, 0.428, 0.996, 0.873, 0.777, 0.945, 0.937, 0.996, 0.59, 0.729, 0.746, 0.925, 0.697, 0.903, 0.841, 0.711, 0.975, 0.654, 0.076, 0.936, 0.525, 0.367, 0.935, 0.74, 0.618, 0.984, 0.971, 0.493, 0.925, 0.831, 0.721, 0.561, 0.69, 0.641, 0.997, 0.908, 0.626, 0.969, 0.907, 0.968, 0.903, 0.973, 0.91, 0.696, 0.597, 0.8, 0.699, 0.405, 0.701, 0.715, 0.584, 0.526, 0.749, 0.959, 0.756, 0.888, 0.602, 0.193, 0.402, 0.29, 0.898, 0.982, 0.691, 0.979, 0.859, 0.65, 0.361, 0.909, 0.447, 0.97, 0.955, 0.993, 0.279, 0.414, 0.66, 0.964, 0.35, 0.999, 0.156, 0.247, 0.589]
global origin = 1
global destination = 50