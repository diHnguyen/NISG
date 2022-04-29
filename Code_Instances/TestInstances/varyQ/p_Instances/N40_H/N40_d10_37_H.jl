global arcs = [1 4; 1 10; 1 20; 1 22; 1 30; 1 32; 2 11; 2 23; 2 33; 3 13; 3 15; 3 24; 3 27; 3 37; 4 13; 4 25; 4 33; 5 6; 5 15; 5 21; 5 26; 6 3; 6 18; 6 24; 6 36; 7 6; 7 25; 8 37; 9 2; 9 10; 9 13; 9 16; 9 22; 9 33; 9 37; 10 16; 10 19; 10 31; 10 37; 11 6; 11 7; 11 25; 11 27; 11 39; 12 6; 12 7; 12 8; 12 23; 12 28; 12 30; 12 34; 13 21; 13 22; 13 30; 14 4; 15 19; 15 25; 15 28; 16 8; 16 40; 17 14; 17 26; 17 31; 17 34; 18 4; 18 5; 18 6; 18 8; 18 22; 18 23; 18 26; 19 33; 20 3; 20 9; 20 10; 20 12; 20 18; 21 4; 21 15; 21 18; 21 26; 21 29; 21 30; 22 2; 22 10; 22 17; 22 27; 23 4; 23 29; 24 22; 25 15; 25 24; 25 39; 26 2; 26 7; 26 19; 26 23; 26 32; 27 17; 27 26; 27 35; 28 5; 28 25; 28 27; 28 36; 29 3; 29 13; 29 30; 29 31; 30 5; 30 17; 30 20; 30 24; 31 18; 32 6; 32 12; 32 13; 32 22; 32 34; 32 37; 33 7; 33 19; 33 24; 33 25; 34 11; 34 18; 34 19; 34 20; 35 5; 35 6; 35 8; 36 3; 36 4; 36 5; 36 16; 36 30; 36 33; 37 3; 37 5; 37 26; 37 38; 38 4; 38 7; 38 18; 38 30; 38 33; 38 36; 38 39; 39 12; 39 17; 39 31; 39 32]
global d_x = [1.0, 3.0, 9.0, 6.0, 10.0, 8.0, 6.0, 8.0, 4.0, 8.0, 9.0, 6.0, 10.0, 9.0, 4.0, 9.0, 6.0, 3.0, 5.0, 10.0, 6.0, 9.0, 7.0, 2.0, 1.0, 1.0, 2.0, 10.0, 3.0, 9.0, 1.0, 6.0, 2.0, 3.0, 1.0, 1.0, 4.0, 8.0, 2.0, 1.0, 7.0, 6.0, 6.0, 6.0, 9.0, 3.0, 3.0, 2.0, 1.0, 3.0, 6.0, 10.0, 3.0, 9.0, 10.0, 3.0, 10.0, 9.0, 3.0, 9.0, 5.0, 6.0, 3.0, 9.0, 7.0, 1.0, 6.0, 4.0, 3.0, 4.0, 5.0, 4.0, 6.0, 4.0, 3.0, 3.0, 3.0, 4.0, 3.0, 5.0, 2.0, 9.0, 7.0, 8.0, 10.0, 3.0, 8.0, 3.0, 9.0, 7.0, 5.0, 7.0, 1.0, 8.0, 3.0, 5.0, 2.0, 8.0, 4.0, 9.0, 1.0, 3.0, 1.0, 1.0, 10.0, 5.0, 1.0, 6.0, 7.0, 9.0, 5.0, 1.0, 2.0, 5.0, 8.0, 1.0, 8.0, 3.0, 8.0, 8.0, 9.0, 3.0, 10.0, 7.0, 8.0, 2.0, 5.0, 7.0, 2.0, 3.0, 2.0, 3.0, 5.0, 3.0, 4.0, 8.0, 2.0, 2.0, 7.0, 5.0, 9.0, 1.0, 2.0, 6.0, 8.0, 10.0, 9.0, 9.0, 10.0, 10.0, 5.0, 7.0]
global b_x = 5
global d_y = [6.0, 2.0, 8.0, 3.0, 7.0, 5.0, 10.0, 1.0, 1.0, 10.0, 7.0, 3.0, 7.0, 7.0, 7.0, 3.0, 9.0, 3.0, 7.0, 8.0, 2.0, 8.0, 5.0, 5.0, 8.0, 5.0, 4.0, 4.0, 9.0, 1.0, 9.0, 10.0, 5.0, 4.0, 9.0, 3.0, 2.0, 4.0, 5.0, 9.0, 10.0, 10.0, 2.0, 6.0, 9.0, 2.0, 9.0, 2.0, 3.0, 8.0, 2.0, 6.0, 5.0, 9.0, 8.0, 6.0, 2.0, 1.0, 10.0, 8.0, 8.0, 2.0, 5.0, 8.0, 6.0, 4.0, 6.0, 4.0, 1.0, 10.0, 2.0, 5.0, 6.0, 5.0, 5.0, 5.0, 1.0, 6.0, 7.0, 4.0, 1.0, 8.0, 4.0, 6.0, 2.0, 2.0, 4.0, 7.0, 10.0, 6.0, 8.0, 4.0, 10.0, 3.0, 2.0, 3.0, 5.0, 6.0, 4.0, 3.0, 6.0, 6.0, 6.0, 7.0, 7.0, 9.0, 2.0, 6.0, 10.0, 5.0, 7.0, 10.0, 8.0, 4.0, 6.0, 2.0, 1.0, 4.0, 4.0, 8.0, 9.0, 5.0, 5.0, 3.0, 9.0, 7.0, 8.0, 9.0, 9.0, 9.0, 10.0, 4.0, 1.0, 8.0, 4.0, 8.0, 6.0, 9.0, 5.0, 7.0, 4.0, 3.0, 4.0, 8.0, 4.0, 4.0, 7.0, 5.0, 8.0, 7.0, 1.0, 4.0]
global b_y = 10
global p = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
global q = [0.812, 0.843, 0.811, 0.922, 0.949, 0.827, 0.844, 0.841, 0.807, 0.969, 0.996, 0.808, 0.926, 0.996, 0.805, 0.935, 0.891, 0.876, 0.912, 0.866, 0.988, 0.995, 0.974, 0.878, 0.974, 0.885, 0.976, 0.805, 0.821, 0.961, 0.944, 0.924, 0.857, 0.859, 0.91, 0.929, 0.823, 0.896, 0.861, 0.951, 0.86, 0.925, 0.934, 0.942, 0.802, 0.92, 0.937, 0.826, 0.92, 0.878, 0.875, 0.814, 0.891, 0.967, 0.933, 0.943, 0.88, 0.995, 0.969, 0.814, 0.926, 0.808, 0.824, 0.915, 0.873, 0.893, 0.97, 0.847, 0.874, 0.851, 0.919, 0.98, 0.985, 0.933, 0.923, 0.867, 0.938, 0.976, 0.947, 0.929, 0.855, 0.994, 0.899, 0.911, 0.847, 0.937, 0.873, 0.803, 0.99, 0.878, 0.918, 0.889, 0.957, 0.872, 0.923, 0.923, 0.825, 0.909, 0.955, 0.918, 0.917, 0.825, 0.976, 0.802, 0.841, 0.82, 0.897, 0.921, 0.848, 0.872, 0.925, 0.975, 0.904, 0.931, 0.806, 0.89, 0.936, 0.843, 0.826, 0.881, 0.977, 0.937, 0.84, 0.915, 0.885, 0.827, 0.925, 0.979, 0.803, 0.942, 0.855, 0.873, 0.828, 0.876, 0.978, 0.866, 0.902, 0.925, 0.925, 0.956, 0.934, 0.803, 0.844, 0.908, 0.867, 0.932, 0.849, 0.819, 0.946, 0.904, 0.961, 0.882]
global origin = 1
global destination = 40