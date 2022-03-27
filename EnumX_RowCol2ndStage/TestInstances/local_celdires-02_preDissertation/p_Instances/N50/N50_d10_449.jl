global arcs = [1 2; 1 4; 1 7; 1 26; 1 31; 1 34; 1 36; 1 46; 2 12; 2 13; 2 18; 2 22; 2 26; 2 46; 2 47; 3 6; 3 12; 3 13; 3 19; 3 34; 3 48; 4 12; 4 32; 4 36; 5 8; 5 15; 6 5; 6 8; 6 10; 6 15; 6 32; 7 16; 7 23; 7 30; 8 30; 8 33; 8 49; 9 2; 9 4; 9 7; 9 13; 9 28; 9 33; 9 35; 9 45; 10 9; 10 33; 10 40; 10 46; 11 21; 11 30; 11 35; 11 47; 12 20; 12 30; 12 34; 13 6; 13 14; 13 19; 13 37; 13 41; 13 49; 14 11; 14 12; 15 14; 15 21; 15 48; 16 7; 16 10; 16 19; 16 27; 16 30; 16 35; 17 22; 17 29; 17 34; 17 36; 17 43; 17 47; 18 10; 18 22; 18 45; 18 46; 19 9; 19 10; 19 24; 20 2; 20 15; 20 22; 20 25; 20 27; 20 28; 21 7; 21 32; 21 43; 22 6; 22 12; 22 17; 22 29; 22 34; 22 35; 22 39; 22 48; 22 50; 23 29; 23 40; 23 43; 24 11; 24 22; 24 50; 25 5; 25 18; 25 30; 25 38; 25 43; 26 4; 26 14; 26 17; 26 43; 27 19; 27 31; 27 37; 28 7; 28 37; 28 45; 29 4; 29 44; 29 50; 30 6; 30 16; 30 26; 30 35; 30 47; 31 3; 31 14; 31 23; 31 28; 31 42; 32 22; 32 28; 32 43; 33 7; 33 14; 33 32; 33 48; 34 5; 34 33; 34 41; 34 43; 34 50; 35 9; 35 11; 35 29; 35 32; 35 46; 35 50; 36 26; 36 37; 36 44; 37 16; 37 24; 37 27; 37 29; 37 47; 38 16; 38 47; 39 17; 39 20; 39 25; 39 40; 39 45; 40 8; 41 7; 41 15; 41 16; 41 40; 41 48; 42 6; 42 19; 42 21; 42 24; 42 49; 43 16; 43 37; 43 44; 44 2; 44 12; 44 24; 45 15; 45 33; 46 6; 46 21; 46 29; 46 31; 47 18; 47 28; 47 44; 48 26; 49 7; 49 36; 49 38; 49 42; 49 47; 49 50]
global d_x = [6.0, 10.0, 3.0, 1.0, 5.0, 4.0, 9.0, 10.0, 7.0, 10.0, 4.0, 2.0, 2.0, 6.0, 9.0, 7.0, 2.0, 9.0, 3.0, 1.0, 3.0, 5.0, 5.0, 10.0, 8.0, 10.0, 5.0, 1.0, 6.0, 9.0, 6.0, 5.0, 7.0, 2.0, 8.0, 6.0, 9.0, 10.0, 9.0, 3.0, 8.0, 6.0, 3.0, 7.0, 10.0, 9.0, 3.0, 9.0, 5.0, 2.0, 9.0, 1.0, 4.0, 5.0, 2.0, 2.0, 6.0, 8.0, 5.0, 3.0, 4.0, 4.0, 2.0, 3.0, 3.0, 9.0, 6.0, 10.0, 4.0, 4.0, 5.0, 9.0, 4.0, 2.0, 10.0, 8.0, 3.0, 1.0, 7.0, 8.0, 10.0, 6.0, 3.0, 7.0, 6.0, 7.0, 5.0, 3.0, 6.0, 4.0, 10.0, 3.0, 6.0, 9.0, 10.0, 9.0, 6.0, 9.0, 4.0, 8.0, 4.0, 2.0, 9.0, 5.0, 9.0, 10.0, 2.0, 1.0, 1.0, 6.0, 6.0, 3.0, 7.0, 6.0, 4.0, 3.0, 10.0, 3.0, 4.0, 6.0, 6.0, 4.0, 9.0, 1.0, 7.0, 5.0, 10.0, 6.0, 9.0, 8.0, 3.0, 5.0, 9.0, 10.0, 6.0, 3.0, 2.0, 8.0, 2.0, 7.0, 3.0, 4.0, 7.0, 1.0, 1.0, 9.0, 4.0, 9.0, 2.0, 6.0, 9.0, 3.0, 8.0, 5.0, 2.0, 5.0, 5.0, 2.0, 3.0, 5.0, 6.0, 6.0, 9.0, 7.0, 5.0, 4.0, 4.0, 10.0, 7.0, 2.0, 9.0, 6.0, 5.0, 7.0, 10.0, 1.0, 7.0, 3.0, 7.0, 5.0, 7.0, 1.0, 3.0, 6.0, 9.0, 7.0, 3.0, 4.0, 6.0, 6.0, 10.0, 6.0, 7.0, 3.0, 2.0, 1.0, 1.0, 8.0, 5.0, 3.0, 10.0, 2.0, 8.0, 8.0]
global b_x = 5
global d_y = [5.0, 3.0, 3.0, 5.0, 3.0, 1.0, 2.0, 6.0, 2.0, 8.0, 1.0, 7.0, 8.0, 2.0, 9.0, 2.0, 9.0, 5.0, 7.0, 5.0, 5.0, 5.0, 2.0, 8.0, 9.0, 9.0, 9.0, 4.0, 3.0, 1.0, 3.0, 3.0, 9.0, 9.0, 7.0, 1.0, 1.0, 9.0, 3.0, 4.0, 10.0, 8.0, 1.0, 4.0, 9.0, 6.0, 5.0, 8.0, 4.0, 5.0, 9.0, 3.0, 2.0, 9.0, 3.0, 1.0, 10.0, 5.0, 7.0, 6.0, 10.0, 5.0, 6.0, 1.0, 5.0, 9.0, 10.0, 4.0, 6.0, 10.0, 3.0, 1.0, 1.0, 2.0, 1.0, 4.0, 4.0, 2.0, 8.0, 5.0, 2.0, 9.0, 9.0, 2.0, 4.0, 10.0, 1.0, 3.0, 4.0, 5.0, 10.0, 6.0, 3.0, 9.0, 7.0, 7.0, 9.0, 4.0, 10.0, 2.0, 1.0, 10.0, 3.0, 7.0, 1.0, 1.0, 3.0, 4.0, 10.0, 4.0, 6.0, 8.0, 8.0, 4.0, 3.0, 5.0, 2.0, 8.0, 10.0, 10.0, 4.0, 8.0, 3.0, 9.0, 7.0, 2.0, 2.0, 7.0, 1.0, 5.0, 6.0, 9.0, 6.0, 4.0, 4.0, 5.0, 3.0, 10.0, 10.0, 5.0, 4.0, 9.0, 5.0, 6.0, 8.0, 5.0, 2.0, 7.0, 10.0, 4.0, 1.0, 1.0, 2.0, 3.0, 5.0, 6.0, 5.0, 2.0, 10.0, 5.0, 6.0, 8.0, 3.0, 10.0, 5.0, 3.0, 9.0, 8.0, 6.0, 10.0, 2.0, 2.0, 4.0, 10.0, 4.0, 8.0, 3.0, 3.0, 9.0, 5.0, 10.0, 4.0, 6.0, 3.0, 8.0, 7.0, 4.0, 8.0, 6.0, 4.0, 8.0, 3.0, 2.0, 4.0, 2.0, 1.0, 7.0, 6.0, 7.0, 4.0, 3.0, 8.0, 7.0, 10.0]
global b_y = 10
global p = [0.187, 0.738, 0.592, 0.846, 0.799, 0.405, 0.467, 0.153, 0.153, 0.491, 0.584, 0.524, 0.004, 0.061, 0.136, 0.679, 0.55, 0.196, 0.702, 0.606, 0.412, 0.646, 0.01, 0.163, 0.104, 0.397, 0.639, 0.046, 0.958, 0.693, 0.93, 0.689, 0.924, 0.034, 0.659, 0.001, 0.165, 0.051, 0.876, 0.647, 0.554, 0.262, 0.195, 0.34, 0.841, 0.422, 0.397, 0.153, 0.613, 0.462, 0.848, 0.797, 0.67, 0.598, 0.362, 0.206, 0.695, 0.559, 0.529, 0.521, 0.351, 0.381, 0.476, 0.449, 0.788, 0.533, 0.678, 0.519, 0.296, 0.559, 0.064, 0.123, 0.188, 0.168, 0.187, 0.312, 0.775, 0.472, 0.142, 0.41, 0.719, 0.743, 0.668, 0.585, 0.266, 0.715, 0.146, 0.678, 0.591, 0.823, 0.864, 0.193, 0.134, 0.146, 0.205, 0.925, 0.592, 0.796, 0.768, 0.598, 0.417, 0.453, 0.703, 0.693, 0.912, 0.719, 0.248, 0.263, 0.201, 0.559, 0.306, 0.918, 0.447, 0.201, 0.999, 0.892, 0.304, 0.054, 0.1, 0.246, 0.926, 0.77, 0.424, 0.705, 0.902, 0.821, 0.71, 0.075, 0.719, 0.547, 0.264, 0.416, 0.653, 0.684, 0.068, 0.231, 0.09, 0.623, 0.841, 0.805, 0.105, 0.214, 0.064, 0.798, 0.521, 0.225, 0.898, 0.809, 0.702, 0.615, 0.228, 0.529, 0.489, 0.165, 0.716, 0.62, 0.635, 0.767, 0.691, 0.04, 0.132, 0.334, 0.033, 0.074, 0.585, 0.961, 0.775, 0.867, 0.744, 0.071, 0.291, 0.613, 0.743, 0.031, 0.3, 0.096, 0.405, 0.103, 0.925, 0.077, 0.504, 0.699, 0.027, 0.211, 0.464, 0.297, 0.01, 0.725, 0.021, 0.687, 0.632, 0.05, 0.272, 0.742, 0.86, 0.355, 0.186, 0.549, 0.44, 0.665, 0.855, 0.711, 0.26, 0.298]
global q = [0.955, 0.769, 0.799, 0.972, 0.947, 0.529, 0.981, 0.283, 0.535, 0.824, 0.605, 0.839, 0.68, 0.978, 0.417, 0.718, 0.557, 0.29, 0.779, 0.763, 0.919, 0.779, 0.618, 0.41, 0.274, 0.528, 0.695, 0.56, 0.96, 0.905, 0.983, 0.724, 0.969, 0.41, 0.922, 0.941, 0.375, 0.359, 0.977, 0.819, 0.766, 0.589, 0.904, 0.759, 0.882, 0.499, 0.886, 0.584, 0.955, 0.665, 0.947, 0.832, 0.971, 0.661, 0.773, 0.289, 0.737, 0.766, 0.696, 0.964, 0.455, 0.879, 0.971, 0.545, 0.793, 0.821, 0.876, 0.672, 0.907, 0.568, 0.943, 0.515, 0.292, 0.734, 0.296, 0.918, 0.955, 0.525, 0.745, 0.511, 0.987, 0.92, 0.779, 0.683, 0.737, 0.768, 0.379, 0.895, 0.626, 0.951, 0.873, 0.382, 0.148, 0.444, 0.757, 0.941, 0.997, 0.944, 0.783, 0.875, 0.987, 0.477, 0.832, 0.737, 0.922, 0.867, 0.683, 0.647, 0.692, 0.824, 0.869, 0.961, 0.679, 0.7, 0.999, 0.931, 0.917, 0.524, 0.382, 0.67, 0.988, 0.874, 0.431, 0.784, 0.971, 0.943, 0.853, 0.667, 0.806, 0.799, 0.475, 0.515, 0.84, 0.827, 0.983, 0.299, 0.819, 0.663, 0.905, 0.974, 0.286, 0.32, 0.595, 0.855, 0.977, 0.371, 0.991, 0.923, 0.806, 0.752, 0.655, 0.784, 0.768, 0.561, 0.728, 0.626, 0.861, 0.778, 0.985, 0.132, 0.235, 0.854, 0.636, 0.777, 0.836, 0.981, 0.862, 0.99, 0.898, 0.357, 0.606, 0.731, 0.995, 0.862, 0.552, 0.91, 0.606, 0.871, 0.987, 0.859, 0.555, 0.838, 0.862, 0.522, 0.569, 0.454, 0.108, 0.819, 0.286, 0.981, 0.882, 0.563, 0.429, 0.862, 0.991, 0.952, 0.335, 0.982, 0.55, 0.94, 0.942, 0.907, 0.477, 0.944]
global origin = 1
global destination = 50