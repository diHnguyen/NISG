global arcs = [1 15; 1 29; 1 32; 1 36; 2 12; 2 15; 2 26; 2 28; 2 29; 3 2; 3 26; 3 43; 4 19; 4 28; 4 32; 4 44; 5 13; 5 18; 5 42; 5 45; 6 7; 6 8; 6 21; 6 22; 6 28; 6 43; 7 3; 7 10; 7 15; 7 16; 7 18; 7 32; 8 7; 8 15; 8 16; 9 13; 9 16; 9 26; 9 32; 9 38; 9 50; 10 13; 10 18; 10 19; 10 39; 10 44; 11 6; 11 19; 11 20; 11 38; 12 6; 12 7; 12 10; 12 23; 13 3; 13 8; 13 18; 13 30; 14 2; 14 17; 14 41; 14 42; 15 3; 15 6; 15 10; 15 19; 15 38; 15 43; 15 46; 16 8; 16 19; 16 22; 16 28; 16 34; 16 46; 16 47; 17 10; 17 25; 17 29; 17 34; 17 41; 17 43; 17 47; 18 4; 18 20; 18 32; 18 38; 18 45; 19 7; 19 11; 19 28; 19 45; 19 46; 20 3; 20 11; 20 19; 20 21; 20 36; 21 14; 21 30; 21 35; 21 44; 21 46; 22 6; 22 37; 23 29; 23 33; 23 35; 23 36; 23 38; 23 46; 24 6; 24 13; 24 16; 24 20; 24 23; 24 26; 24 29; 25 27; 25 35; 26 4; 26 11; 26 32; 26 33; 26 42; 26 45; 27 10; 27 38; 27 49; 28 18; 28 21; 28 44; 28 46; 29 11; 29 20; 29 22; 29 40; 29 42; 29 46; 30 14; 30 18; 30 43; 30 46; 31 7; 31 10; 31 16; 31 21; 31 40; 32 16; 32 21; 32 33; 32 50; 33 21; 33 22; 33 32; 34 50; 35 9; 35 18; 35 19; 36 27; 36 41; 37 5; 37 20; 37 26; 37 41; 37 50; 38 32; 38 39; 38 46; 38 47; 39 8; 39 21; 39 29; 39 32; 39 42; 39 45; 40 31; 40 37; 40 39; 40 41; 40 43; 40 44; 41 3; 41 18; 41 22; 41 24; 41 50; 42 3; 42 25; 42 36; 42 46; 43 14; 43 24; 43 30; 43 41; 44 2; 44 4; 44 7; 44 25; 44 47; 45 4; 45 16; 45 24; 45 30; 45 33; 45 38; 45 42; 46 13; 46 17; 46 23; 46 36; 47 4; 47 9; 47 23; 47 24; 47 48; 48 5; 48 22; 48 34; 49 6; 49 7; 49 9; 49 12; 49 26; 49 39]
global d_x = [8.0, 2.0, 9.0, 1.0, 6.0, 5.0, 5.0, 6.0, 9.0, 4.0, 6.0, 7.0, 9.0, 3.0, 10.0, 8.0, 7.0, 2.0, 7.0, 2.0, 3.0, 8.0, 5.0, 10.0, 2.0, 10.0, 7.0, 4.0, 6.0, 3.0, 7.0, 4.0, 2.0, 1.0, 9.0, 8.0, 7.0, 6.0, 3.0, 6.0, 9.0, 9.0, 10.0, 9.0, 8.0, 8.0, 4.0, 8.0, 7.0, 9.0, 6.0, 7.0, 7.0, 7.0, 4.0, 3.0, 7.0, 8.0, 9.0, 5.0, 10.0, 4.0, 5.0, 8.0, 4.0, 7.0, 10.0, 4.0, 10.0, 10.0, 9.0, 9.0, 9.0, 2.0, 6.0, 8.0, 6.0, 4.0, 6.0, 6.0, 6.0, 6.0, 2.0, 7.0, 10.0, 6.0, 8.0, 5.0, 2.0, 2.0, 2.0, 2.0, 7.0, 1.0, 3.0, 5.0, 10.0, 4.0, 8.0, 8.0, 1.0, 6.0, 4.0, 3.0, 8.0, 7.0, 9.0, 3.0, 8.0, 8.0, 9.0, 9.0, 5.0, 2.0, 3.0, 10.0, 5.0, 6.0, 1.0, 8.0, 4.0, 8.0, 10.0, 10.0, 1.0, 7.0, 9.0, 8.0, 7.0, 5.0, 10.0, 6.0, 9.0, 5.0, 8.0, 6.0, 3.0, 10.0, 2.0, 3.0, 8.0, 7.0, 10.0, 10.0, 5.0, 9.0, 9.0, 1.0, 3.0, 1.0, 5.0, 6.0, 3.0, 3.0, 10.0, 3.0, 1.0, 10.0, 2.0, 9.0, 6.0, 5.0, 4.0, 3.0, 1.0, 10.0, 10.0, 1.0, 6.0, 10.0, 4.0, 10.0, 8.0, 6.0, 8.0, 9.0, 4.0, 3.0, 3.0, 6.0, 8.0, 1.0, 7.0, 7.0, 3.0, 7.0, 10.0, 2.0, 4.0, 10.0, 4.0, 1.0, 1.0, 7.0, 3.0, 5.0, 9.0, 6.0, 1.0, 7.0, 6.0, 1.0, 5.0, 5.0, 9.0, 2.0, 9.0, 9.0, 5.0, 9.0, 7.0, 10.0, 3.0, 4.0, 10.0, 3.0, 3.0, 1.0, 7.0, 2.0, 1.0, 5.0, 1.0, 8.0, 3.0]
global b_x = 5
global d_y = [4.0, 10.0, 1.0, 6.0, 10.0, 8.0, 7.0, 5.0, 8.0, 8.0, 2.0, 9.0, 2.0, 7.0, 1.0, 1.0, 2.0, 6.0, 1.0, 10.0, 2.0, 4.0, 2.0, 5.0, 1.0, 10.0, 7.0, 2.0, 9.0, 5.0, 7.0, 1.0, 7.0, 1.0, 7.0, 7.0, 10.0, 5.0, 5.0, 9.0, 5.0, 5.0, 1.0, 8.0, 3.0, 4.0, 10.0, 2.0, 6.0, 8.0, 3.0, 3.0, 2.0, 1.0, 6.0, 6.0, 1.0, 9.0, 6.0, 1.0, 8.0, 10.0, 7.0, 2.0, 6.0, 4.0, 2.0, 4.0, 4.0, 3.0, 10.0, 1.0, 6.0, 6.0, 10.0, 2.0, 1.0, 7.0, 3.0, 6.0, 2.0, 9.0, 7.0, 10.0, 2.0, 2.0, 3.0, 10.0, 10.0, 2.0, 7.0, 4.0, 2.0, 5.0, 4.0, 6.0, 5.0, 9.0, 5.0, 4.0, 9.0, 5.0, 9.0, 4.0, 6.0, 3.0, 10.0, 3.0, 8.0, 6.0, 7.0, 3.0, 5.0, 9.0, 6.0, 3.0, 9.0, 3.0, 10.0, 4.0, 7.0, 8.0, 7.0, 1.0, 1.0, 5.0, 2.0, 1.0, 7.0, 8.0, 7.0, 5.0, 1.0, 8.0, 5.0, 3.0, 8.0, 6.0, 1.0, 6.0, 10.0, 9.0, 8.0, 5.0, 3.0, 7.0, 7.0, 7.0, 3.0, 4.0, 4.0, 2.0, 5.0, 5.0, 4.0, 3.0, 2.0, 4.0, 2.0, 10.0, 5.0, 9.0, 3.0, 7.0, 3.0, 3.0, 5.0, 3.0, 10.0, 9.0, 3.0, 5.0, 8.0, 5.0, 5.0, 6.0, 3.0, 5.0, 9.0, 3.0, 7.0, 6.0, 6.0, 9.0, 10.0, 7.0, 3.0, 4.0, 4.0, 2.0, 4.0, 10.0, 1.0, 2.0, 4.0, 2.0, 9.0, 7.0, 5.0, 7.0, 9.0, 2.0, 10.0, 8.0, 5.0, 9.0, 3.0, 9.0, 1.0, 9.0, 10.0, 3.0, 7.0, 9.0, 1.0, 10.0, 2.0, 5.0, 2.0, 6.0, 4.0, 9.0, 6.0, 2.0, 2.0]
global b_y = 10
global p = [0.645, 0.945, 0.014, 0.598, 0.1, 0.577, 0.055, 0.349, 0.653, 0.274, 0.311, 0.374, 0.503, 0.975, 0.684, 0.741, 0.25, 0.513, 0.389, 0.121, 0.085, 0.168, 0.158, 0.665, 0.806, 0.755, 0.485, 0.887, 0.453, 0.908, 0.428, 0.352, 0.48, 0.653, 0.829, 0.845, 0.834, 0.199, 0.308, 0.162, 0.099, 0.657, 0.895, 0.453, 0.884, 0.967, 0.048, 0.479, 0.489, 0.576, 0.538, 0.345, 0.622, 0.453, 0.128, 0.26, 0.208, 0.179, 0.851, 0.926, 0.891, 0.312, 0.199, 0.915, 0.081, 0.205, 0.137, 0.197, 0.388, 0.802, 0.528, 0.386, 0.553, 0.263, 0.968, 0.819, 0.038, 0.57, 0.433, 0.662, 0.785, 0.693, 0.495, 0.625, 0.261, 0.964, 0.585, 0.38, 0.62, 0.837, 0.087, 0.87, 0.284, 0.514, 0.367, 0.919, 0.542, 0.312, 0.484, 0.606, 0.999, 0.632, 0.615, 0.512, 0.2, 0.756, 0.615, 0.715, 0.279, 0.455, 0.368, 0.694, 0.388, 0.175, 0.066, 0.092, 0.407, 0.11, 0.816, 0.582, 0.677, 0.143, 0.158, 0.454, 0.396, 0.246, 0.822, 0.982, 0.939, 0.741, 0.812, 0.391, 0.255, 0.628, 0.109, 0.018, 0.009, 0.767, 0.789, 0.235, 0.187, 0.004, 0.685, 0.174, 0.46, 0.681, 0.422, 0.522, 0.762, 0.185, 0.071, 0.539, 0.138, 0.417, 0.103, 0.388, 0.436, 0.165, 0.781, 0.141, 0.298, 0.657, 0.28, 0.677, 0.106, 0.271, 0.966, 0.051, 0.577, 0.211, 0.065, 0.063, 0.35, 0.977, 0.283, 0.581, 0.567, 0.434, 0.5, 0.355, 0.084, 0.202, 0.535, 0.807, 0.336, 0.094, 0.175, 0.006, 0.173, 0.58, 0.951, 0.638, 0.431, 0.902, 0.973, 0.702, 0.946, 0.184, 0.178, 0.859, 0.392, 0.966, 0.359, 0.75, 0.652, 0.631, 0.171, 0.04, 0.947, 0.279, 0.674, 0.418, 0.048, 0.637, 0.581, 0.818, 0.169, 0.015, 0.664, 0.892, 0.039, 0.804, 0.655, 0.294, 0.783]
global q = [0.807, 0.995, 0.042, 0.723, 0.677, 0.919, 0.933, 0.772, 0.705, 0.58, 0.646, 0.549, 0.548, 0.996, 0.95, 0.865, 0.331, 0.643, 0.647, 0.679, 0.333, 0.257, 0.824, 0.683, 0.869, 0.827, 0.733, 0.896, 0.77, 0.968, 0.724, 0.526, 0.876, 0.955, 0.908, 0.895, 0.95, 0.578, 0.364, 0.459, 0.224, 0.971, 0.955, 0.671, 0.93, 0.982, 0.626, 0.828, 0.736, 0.985, 0.65, 0.54, 0.972, 0.679, 0.763, 0.397, 0.721, 0.789, 0.992, 0.972, 0.983, 0.963, 0.578, 0.991, 0.176, 0.253, 0.632, 0.832, 0.511, 0.862, 0.792, 0.479, 0.779, 0.293, 0.97, 0.925, 0.563, 0.746, 0.814, 0.832, 0.838, 0.837, 0.815, 0.683, 0.821, 0.99, 0.885, 0.526, 0.64, 0.913, 0.231, 0.9, 0.537, 0.536, 0.838, 0.972, 0.606, 0.407, 0.903, 0.737, 0.999, 0.798, 0.714, 0.72, 0.427, 0.984, 0.663, 0.733, 0.847, 0.765, 0.46, 0.984, 0.819, 0.872, 0.472, 0.487, 0.502, 0.99, 0.955, 0.793, 0.908, 0.183, 0.292, 0.739, 0.753, 0.623, 0.894, 0.996, 0.97, 0.816, 0.993, 0.65, 0.981, 0.913, 0.235, 0.926, 0.124, 0.941, 0.879, 0.889, 0.443, 0.513, 0.996, 0.73, 0.855, 0.895, 0.743, 0.981, 0.769, 0.917, 0.677, 0.641, 0.304, 0.443, 0.463, 0.901, 0.775, 0.82, 0.871, 0.976, 0.536, 0.996, 0.443, 0.721, 0.492, 0.37, 0.98, 0.429, 0.91, 0.817, 0.99, 0.856, 0.91, 0.994, 0.734, 0.828, 0.808, 0.841, 0.849, 0.561, 0.391, 0.243, 0.854, 0.979, 0.999, 0.562, 0.667, 0.187, 0.418, 0.74, 0.976, 0.714, 0.735, 0.974, 0.975, 0.729, 0.971, 0.622, 0.99, 0.901, 0.759, 0.966, 0.547, 0.896, 0.849, 0.851, 0.295, 0.689, 0.965, 0.726, 0.756, 0.862, 0.668, 0.671, 0.767, 0.985, 0.82, 0.505, 0.76, 0.896, 0.985, 0.964, 0.819, 0.357, 0.836]
global origin = 1
global destination = 50