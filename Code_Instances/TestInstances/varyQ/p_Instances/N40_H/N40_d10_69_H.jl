global arcs = [1 15; 1 19; 1 33; 1 37; 2 6; 2 39; 3 12; 3 25; 3 38; 4 10; 4 11; 4 18; 4 26; 4 38; 4 39; 5 9; 5 13; 5 21; 5 22; 5 26; 5 40; 6 11; 6 17; 6 26; 7 2; 7 10; 7 15; 7 21; 7 32; 7 35; 8 6; 8 19; 8 21; 8 22; 8 31; 8 36; 9 15; 9 22; 9 30; 10 19; 10 27; 10 32; 11 10; 11 12; 11 20; 11 35; 12 10; 12 18; 12 23; 12 28; 13 19; 13 22; 13 29; 14 12; 14 25; 14 35; 14 40; 15 34; 15 39; 16 6; 16 22; 17 21; 17 32; 18 2; 18 27; 18 28; 18 33; 19 7; 19 15; 19 16; 19 17; 19 28; 19 29; 19 32; 20 4; 20 30; 21 4; 21 9; 21 13; 21 17; 21 34; 22 2; 22 3; 22 4; 22 5; 22 6; 22 13; 22 18; 22 34; 22 40; 23 11; 23 16; 23 33; 23 36; 23 40; 24 2; 24 13; 24 16; 25 5; 25 13; 25 26; 25 34; 25 35; 25 38; 26 2; 26 23; 26 25; 26 31; 27 8; 27 21; 28 10; 28 11; 28 17; 28 31; 28 39; 29 9; 30 8; 30 26; 31 9; 31 25; 31 26; 31 32; 31 34; 32 3; 32 24; 32 39; 33 18; 33 36; 33 39; 34 4; 34 9; 34 18; 34 19; 34 31; 35 5; 35 7; 35 12; 35 14; 35 15; 35 19; 35 38; 36 11; 36 12; 36 25; 36 26; 36 34; 37 7; 37 9; 37 19; 37 23; 38 2; 38 20; 39 10; 39 11; 39 19; 39 24]
global d_x = [5.0, 7.0, 1.0, 3.0, 9.0, 2.0, 10.0, 5.0, 5.0, 2.0, 8.0, 7.0, 9.0, 8.0, 10.0, 6.0, 5.0, 8.0, 5.0, 2.0, 3.0, 8.0, 8.0, 10.0, 2.0, 4.0, 7.0, 9.0, 4.0, 5.0, 6.0, 4.0, 10.0, 2.0, 9.0, 2.0, 5.0, 8.0, 7.0, 8.0, 2.0, 9.0, 8.0, 3.0, 9.0, 6.0, 4.0, 2.0, 4.0, 9.0, 7.0, 8.0, 5.0, 5.0, 3.0, 10.0, 4.0, 5.0, 8.0, 3.0, 7.0, 3.0, 5.0, 1.0, 6.0, 8.0, 4.0, 2.0, 8.0, 7.0, 8.0, 10.0, 3.0, 8.0, 2.0, 8.0, 7.0, 9.0, 2.0, 4.0, 4.0, 6.0, 7.0, 5.0, 3.0, 2.0, 8.0, 1.0, 5.0, 2.0, 5.0, 7.0, 7.0, 7.0, 2.0, 1.0, 7.0, 6.0, 10.0, 9.0, 7.0, 9.0, 3.0, 1.0, 1.0, 5.0, 8.0, 1.0, 3.0, 4.0, 8.0, 3.0, 6.0, 5.0, 1.0, 4.0, 6.0, 9.0, 2.0, 6.0, 1.0, 10.0, 5.0, 6.0, 1.0, 8.0, 8.0, 10.0, 1.0, 5.0, 6.0, 4.0, 1.0, 8.0, 8.0, 1.0, 6.0, 1.0, 9.0, 5.0, 5.0, 8.0, 6.0, 6.0, 10.0, 9.0, 5.0, 2.0, 6.0, 5.0, 4.0, 1.0, 9.0, 6.0, 8.0, 7.0]
global b_x = 5
global d_y = [3.0, 1.0, 9.0, 9.0, 3.0, 8.0, 4.0, 10.0, 8.0, 9.0, 8.0, 9.0, 6.0, 7.0, 5.0, 7.0, 7.0, 4.0, 9.0, 9.0, 2.0, 3.0, 2.0, 3.0, 5.0, 9.0, 4.0, 9.0, 8.0, 2.0, 9.0, 1.0, 4.0, 1.0, 7.0, 6.0, 4.0, 3.0, 7.0, 4.0, 4.0, 10.0, 1.0, 3.0, 10.0, 9.0, 4.0, 1.0, 9.0, 3.0, 10.0, 10.0, 1.0, 6.0, 10.0, 1.0, 2.0, 1.0, 2.0, 6.0, 2.0, 3.0, 8.0, 4.0, 7.0, 2.0, 5.0, 3.0, 1.0, 3.0, 4.0, 8.0, 10.0, 3.0, 1.0, 1.0, 3.0, 1.0, 6.0, 10.0, 6.0, 4.0, 8.0, 6.0, 8.0, 4.0, 7.0, 10.0, 5.0, 1.0, 2.0, 9.0, 4.0, 3.0, 6.0, 9.0, 6.0, 4.0, 3.0, 3.0, 8.0, 6.0, 2.0, 3.0, 1.0, 8.0, 6.0, 9.0, 3.0, 1.0, 7.0, 7.0, 9.0, 2.0, 6.0, 6.0, 3.0, 4.0, 2.0, 6.0, 6.0, 10.0, 6.0, 5.0, 3.0, 3.0, 2.0, 4.0, 10.0, 8.0, 2.0, 7.0, 9.0, 5.0, 8.0, 4.0, 1.0, 4.0, 4.0, 2.0, 10.0, 8.0, 8.0, 2.0, 10.0, 5.0, 1.0, 5.0, 10.0, 8.0, 6.0, 7.0, 5.0, 8.0, 2.0, 8.0]
global b_y = 10
global p = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
global q = [0.823, 0.925, 0.994, 0.983, 0.835, 0.976, 0.964, 0.813, 0.968, 0.989, 0.949, 0.859, 0.914, 0.905, 0.911, 0.952, 0.848, 0.952, 0.96, 0.856, 0.962, 0.928, 0.969, 0.985, 0.926, 0.883, 0.863, 0.984, 0.991, 0.916, 0.847, 0.897, 0.927, 0.944, 0.821, 0.826, 0.814, 0.804, 0.969, 0.932, 0.876, 0.924, 0.912, 0.883, 0.943, 0.984, 0.811, 0.892, 0.981, 0.836, 0.971, 0.96, 0.988, 0.812, 0.909, 0.897, 0.869, 0.801, 0.98, 0.909, 0.995, 0.801, 0.903, 0.838, 0.969, 0.96, 0.989, 0.901, 0.989, 0.961, 0.815, 0.898, 0.896, 0.891, 0.922, 0.861, 0.804, 0.869, 0.812, 0.962, 0.819, 0.865, 0.971, 0.8, 0.96, 0.973, 0.932, 0.938, 0.868, 0.914, 0.86, 0.924, 0.869, 0.94, 0.834, 0.806, 0.904, 0.966, 0.934, 0.909, 0.805, 0.975, 0.959, 0.988, 0.964, 0.88, 0.912, 0.863, 0.905, 0.819, 0.938, 0.859, 0.902, 0.898, 0.807, 0.853, 0.822, 0.922, 0.919, 0.967, 0.86, 0.951, 0.995, 0.89, 0.841, 0.989, 0.812, 0.838, 0.967, 0.996, 0.95, 0.942, 0.881, 0.923, 0.981, 0.872, 0.858, 0.988, 0.931, 0.972, 0.947, 0.816, 0.955, 0.938, 0.819, 0.955, 0.892, 0.956, 0.963, 0.883, 0.835, 0.983, 0.892, 0.989, 0.84, 0.873]
global origin = 1
global destination = 40