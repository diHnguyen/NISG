global arcs = [1 2; 1 14; 1 21; 1 26; 1 27; 1 31; 1 33; 1 37; 2 8; 3 11; 4 5; 4 12; 4 22; 5 8; 5 13; 5 37; 6 5; 6 16; 6 25; 7 5; 7 13; 7 15; 7 34; 8 11; 8 16; 8 40; 9 4; 9 17; 9 21; 9 22; 9 29; 9 33; 9 37; 9 38; 10 5; 10 8; 10 15; 10 19; 10 20; 10 33; 10 40; 11 4; 11 5; 11 12; 11 13; 11 24; 11 28; 11 36; 11 37; 12 4; 12 6; 12 19; 12 22; 13 10; 13 24; 13 25; 13 34; 14 9; 14 15; 14 39; 15 4; 15 7; 15 25; 15 28; 15 36; 16 13; 16 18; 16 20; 16 30; 16 36; 17 12; 17 27; 18 23; 19 6; 19 33; 19 38; 20 12; 20 13; 20 21; 20 28; 20 29; 20 35; 21 7; 21 13; 22 6; 22 7; 22 8; 22 18; 22 21; 22 23; 22 31; 23 5; 23 7; 23 14; 23 29; 23 30; 23 31; 24 8; 24 13; 24 22; 24 39; 25 7; 25 8; 25 17; 25 22; 25 23; 25 35; 25 38; 26 22; 26 23; 26 25; 26 33; 26 36; 27 18; 27 23; 28 3; 28 21; 28 22; 28 29; 28 35; 29 15; 29 16; 29 32; 30 6; 30 8; 30 9; 30 35; 30 37; 30 40; 31 2; 31 13; 32 4; 32 10; 32 14; 32 19; 32 27; 32 28; 32 33; 32 38; 33 4; 33 5; 33 8; 33 14; 33 21; 34 2; 34 3; 34 13; 34 22; 34 40; 35 2; 35 30; 35 32; 35 38; 36 7; 36 10; 36 19; 36 27; 36 28; 36 29; 36 31; 36 35; 37 2; 37 11; 37 31; 38 3; 38 10; 39 4; 39 29; 39 35]
global d_x = [4.0, 4.0, 1.0, 1.0, 5.0, 1.0, 1.0, 8.0, 4.0, 10.0, 6.0, 7.0, 3.0, 3.0, 9.0, 3.0, 1.0, 9.0, 8.0, 2.0, 10.0, 1.0, 6.0, 4.0, 4.0, 2.0, 6.0, 5.0, 1.0, 3.0, 5.0, 6.0, 8.0, 3.0, 8.0, 6.0, 3.0, 9.0, 7.0, 6.0, 8.0, 5.0, 7.0, 9.0, 10.0, 5.0, 4.0, 10.0, 7.0, 1.0, 5.0, 8.0, 8.0, 5.0, 7.0, 6.0, 4.0, 4.0, 9.0, 9.0, 4.0, 3.0, 2.0, 2.0, 3.0, 3.0, 10.0, 7.0, 7.0, 3.0, 1.0, 6.0, 9.0, 7.0, 7.0, 9.0, 10.0, 8.0, 9.0, 7.0, 5.0, 9.0, 10.0, 4.0, 1.0, 6.0, 10.0, 3.0, 8.0, 5.0, 5.0, 9.0, 3.0, 10.0, 1.0, 4.0, 4.0, 1.0, 3.0, 10.0, 8.0, 1.0, 1.0, 7.0, 9.0, 9.0, 2.0, 5.0, 5.0, 5.0, 7.0, 9.0, 3.0, 10.0, 5.0, 2.0, 7.0, 6.0, 10.0, 7.0, 8.0, 9.0, 4.0, 8.0, 7.0, 10.0, 3.0, 10.0, 2.0, 2.0, 2.0, 8.0, 4.0, 8.0, 10.0, 6.0, 2.0, 7.0, 10.0, 1.0, 2.0, 6.0, 9.0, 3.0, 9.0, 4.0, 10.0, 4.0, 10.0, 5.0, 9.0, 4.0, 7.0, 8.0, 8.0, 1.0, 8.0, 7.0, 5.0, 9.0, 1.0, 5.0, 4.0, 5.0, 1.0, 3.0, 10.0, 9.0, 2.0]
global b_x = 5
global d_y = [4.0, 6.0, 7.0, 9.0, 2.0, 6.0, 7.0, 1.0, 10.0, 6.0, 9.0, 8.0, 2.0, 4.0, 2.0, 7.0, 9.0, 1.0, 5.0, 5.0, 3.0, 9.0, 3.0, 6.0, 10.0, 2.0, 7.0, 5.0, 9.0, 3.0, 2.0, 6.0, 10.0, 2.0, 2.0, 7.0, 1.0, 5.0, 1.0, 9.0, 9.0, 6.0, 3.0, 1.0, 9.0, 10.0, 5.0, 2.0, 7.0, 5.0, 6.0, 4.0, 1.0, 3.0, 5.0, 7.0, 9.0, 4.0, 4.0, 4.0, 7.0, 3.0, 1.0, 1.0, 8.0, 5.0, 5.0, 8.0, 10.0, 5.0, 4.0, 3.0, 1.0, 1.0, 10.0, 4.0, 7.0, 1.0, 1.0, 5.0, 3.0, 2.0, 6.0, 10.0, 10.0, 7.0, 9.0, 9.0, 1.0, 10.0, 2.0, 6.0, 8.0, 7.0, 8.0, 5.0, 10.0, 2.0, 6.0, 5.0, 6.0, 2.0, 2.0, 9.0, 5.0, 5.0, 1.0, 9.0, 10.0, 10.0, 9.0, 5.0, 3.0, 9.0, 9.0, 6.0, 4.0, 2.0, 1.0, 5.0, 10.0, 8.0, 2.0, 8.0, 3.0, 3.0, 6.0, 5.0, 4.0, 8.0, 10.0, 10.0, 8.0, 6.0, 10.0, 6.0, 1.0, 3.0, 4.0, 3.0, 4.0, 4.0, 8.0, 2.0, 9.0, 4.0, 10.0, 1.0, 9.0, 8.0, 1.0, 3.0, 6.0, 1.0, 8.0, 7.0, 6.0, 1.0, 1.0, 3.0, 10.0, 7.0, 9.0, 2.0, 4.0, 10.0, 7.0, 4.0, 2.0]
global b_y = 10
global p = [0.441, 0.211, 0.964, 0.228, 0.748, 0.686, 0.45, 0.775, 0.095, 0.883, 0.639, 0.539, 0.518, 0.487, 0.135, 0.774, 0.631, 0.421, 0.615, 0.296, 0.524, 0.35, 0.581, 0.681, 0.363, 0.475, 0.398, 0.658, 0.429, 0.546, 0.475, 0.597, 0.572, 0.385, 0.759, 0.42, 0.717, 0.567, 0.949, 0.035, 0.072, 0.274, 0.271, 0.579, 0.014, 0.339, 0.695, 0.475, 0.691, 0.073, 0.489, 0.962, 0.242, 0.055, 0.055, 0.515, 0.309, 0.561, 0.153, 0.116, 0.856, 0.719, 0.941, 0.608, 0.559, 0.209, 0.024, 0.021, 0.678, 0.155, 0.524, 0.068, 0.963, 0.67, 0.628, 0.752, 0.571, 0.201, 0.24, 0.894, 0.437, 0.393, 0.315, 0.463, 0.486, 0.391, 0.571, 0.234, 0.294, 0.746, 0.718, 0.002, 0.018, 0.543, 0.988, 0.052, 0.384, 0.002, 0.764, 0.487, 0.303, 0.665, 0.99, 0.051, 0.757, 0.48, 0.23, 0.198, 0.517, 0.915, 0.627, 0.945, 0.976, 0.181, 0.09, 0.661, 0.294, 0.709, 0.741, 0.569, 0.9, 0.186, 0.085, 0.778, 0.69, 0.425, 0.863, 0.466, 0.773, 0.889, 0.42, 0.263, 0.38, 0.605, 0.265, 0.989, 0.212, 0.726, 0.033, 0.29, 0.098, 0.983, 0.175, 0.271, 0.661, 0.483, 0.38, 0.224, 0.022, 0.255, 0.16, 0.012, 0.289, 0.167, 0.47, 0.59, 0.286, 0.526, 0.377, 0.239, 0.757, 0.139, 0.214, 0.941, 0.482, 0.697, 0.437, 0.947, 0.508]
global q = [0.599, 0.785, 0.97, 0.932, 0.887, 0.848, 0.456, 0.962, 0.978, 0.9, 0.774, 0.587, 0.925, 0.653, 0.813, 0.782, 0.865, 0.518, 0.925, 0.857, 0.875, 0.644, 0.917, 0.912, 0.645, 0.734, 0.88, 0.693, 0.474, 0.718, 0.756, 0.731, 0.881, 0.776, 0.987, 0.692, 0.764, 0.728, 0.964, 0.138, 0.545, 0.952, 0.968, 0.774, 0.456, 0.557, 0.711, 0.936, 0.81, 0.938, 0.884, 0.987, 0.903, 0.443, 0.726, 0.671, 0.319, 0.903, 0.741, 0.991, 0.942, 0.825, 0.992, 0.791, 0.77, 0.646, 0.065, 0.266, 0.823, 0.967, 0.694, 0.806, 0.98, 0.774, 0.889, 0.858, 0.777, 0.287, 0.876, 0.929, 0.624, 0.491, 0.394, 0.677, 0.978, 0.818, 0.731, 0.741, 0.604, 0.834, 0.824, 0.437, 0.308, 0.985, 0.993, 0.44, 0.955, 0.908, 0.918, 0.579, 0.896, 0.698, 0.997, 0.979, 0.82, 0.747, 0.267, 0.913, 0.562, 0.945, 0.997, 0.999, 0.981, 0.769, 0.489, 0.979, 0.718, 0.863, 0.855, 0.581, 0.974, 0.842, 0.576, 0.793, 0.983, 0.643, 0.973, 0.566, 0.901, 0.947, 0.435, 0.741, 0.456, 0.763, 0.615, 0.999, 0.82, 0.923, 0.67, 0.359, 0.605, 0.986, 0.76, 0.523, 0.826, 0.8, 0.414, 0.289, 0.049, 0.334, 0.715, 0.477, 0.705, 0.419, 0.721, 0.592, 0.349, 0.909, 0.844, 0.796, 0.768, 0.923, 0.589, 0.965, 0.943, 0.986, 0.723, 0.962, 0.547]
global origin = 1
global destination = 40