global arcs = [1 10; 1 15; 1 20; 1 25; 1 37; 2 5; 2 8; 2 15; 2 36; 3 5; 3 9; 3 12; 3 17; 3 20; 3 23; 3 26; 4 9; 4 10; 4 11; 4 30; 4 36; 5 6; 5 28; 6 4; 6 28; 7 11; 7 12; 7 30; 7 38; 8 10; 8 13; 8 31; 8 38; 9 5; 9 22; 9 35; 9 36; 10 34; 10 36; 11 15; 11 24; 11 31; 11 32; 11 39; 12 4; 12 17; 12 23; 12 33; 12 39; 13 8; 13 19; 13 31; 14 25; 14 26; 14 31; 14 34; 15 3; 15 23; 15 28; 16 15; 16 29; 17 3; 17 9; 17 15; 18 4; 18 7; 18 17; 18 19; 18 30; 18 34; 18 38; 19 2; 19 21; 19 33; 19 37; 20 10; 20 15; 20 31; 21 9; 21 13; 22 21; 23 27; 23 34; 24 3; 24 11; 24 27; 25 14; 25 20; 25 21; 25 23; 26 4; 26 8; 26 22; 26 31; 26 32; 27 9; 27 21; 28 5; 28 9; 28 20; 28 26; 28 34; 28 36; 29 3; 29 10; 29 18; 29 19; 29 20; 29 21; 29 27; 29 32; 29 37; 29 40; 30 13; 30 22; 30 26; 30 28; 30 36; 31 5; 31 19; 31 33; 31 36; 32 5; 32 8; 32 14; 32 38; 33 6; 33 7; 33 9; 33 19; 33 23; 33 24; 34 11; 34 18; 34 29; 34 33; 34 40; 35 3; 35 10; 35 14; 35 17; 35 27; 35 28; 36 13; 36 21; 36 22; 36 29; 36 33; 37 13; 37 16; 38 3; 38 6; 38 26; 38 40; 39 12; 39 14; 39 28; 39 34]
global d_x = [10.0, 3.0, 9.0, 3.0, 2.0, 7.0, 1.0, 1.0, 10.0, 7.0, 5.0, 9.0, 7.0, 5.0, 8.0, 9.0, 5.0, 6.0, 7.0, 8.0, 7.0, 1.0, 6.0, 5.0, 6.0, 10.0, 3.0, 10.0, 3.0, 5.0, 10.0, 2.0, 5.0, 9.0, 9.0, 4.0, 4.0, 9.0, 10.0, 1.0, 1.0, 3.0, 4.0, 2.0, 7.0, 7.0, 4.0, 3.0, 9.0, 9.0, 2.0, 9.0, 2.0, 3.0, 8.0, 8.0, 2.0, 6.0, 1.0, 2.0, 4.0, 3.0, 8.0, 10.0, 5.0, 7.0, 6.0, 4.0, 8.0, 10.0, 4.0, 7.0, 4.0, 10.0, 9.0, 10.0, 9.0, 6.0, 9.0, 10.0, 6.0, 2.0, 2.0, 8.0, 6.0, 4.0, 10.0, 6.0, 3.0, 8.0, 7.0, 7.0, 1.0, 1.0, 10.0, 8.0, 5.0, 10.0, 8.0, 4.0, 4.0, 4.0, 1.0, 6.0, 8.0, 5.0, 8.0, 2.0, 8.0, 10.0, 10.0, 10.0, 3.0, 2.0, 6.0, 7.0, 3.0, 8.0, 4.0, 7.0, 8.0, 3.0, 1.0, 6.0, 4.0, 3.0, 10.0, 6.0, 5.0, 6.0, 5.0, 10.0, 2.0, 1.0, 4.0, 10.0, 1.0, 2.0, 8.0, 7.0, 10.0, 9.0, 5.0, 6.0, 2.0, 2.0, 10.0, 9.0, 9.0, 10.0, 9.0, 9.0, 1.0, 6.0, 8.0, 9.0, 1.0, 2.0]
global b_x = 5
global d_y = [9.0, 10.0, 9.0, 1.0, 2.0, 3.0, 4.0, 6.0, 2.0, 9.0, 3.0, 9.0, 10.0, 6.0, 10.0, 3.0, 10.0, 6.0, 1.0, 4.0, 7.0, 5.0, 8.0, 1.0, 3.0, 10.0, 9.0, 7.0, 3.0, 9.0, 1.0, 1.0, 1.0, 1.0, 5.0, 2.0, 5.0, 1.0, 5.0, 10.0, 5.0, 10.0, 4.0, 4.0, 3.0, 2.0, 5.0, 8.0, 4.0, 5.0, 7.0, 10.0, 6.0, 1.0, 10.0, 5.0, 2.0, 6.0, 4.0, 3.0, 3.0, 7.0, 7.0, 3.0, 6.0, 4.0, 1.0, 1.0, 7.0, 1.0, 1.0, 8.0, 3.0, 8.0, 2.0, 5.0, 8.0, 7.0, 3.0, 8.0, 3.0, 3.0, 1.0, 6.0, 2.0, 5.0, 9.0, 7.0, 2.0, 2.0, 1.0, 3.0, 3.0, 5.0, 8.0, 4.0, 1.0, 8.0, 7.0, 6.0, 10.0, 8.0, 2.0, 3.0, 6.0, 6.0, 6.0, 8.0, 6.0, 7.0, 10.0, 9.0, 3.0, 1.0, 7.0, 2.0, 2.0, 8.0, 8.0, 5.0, 9.0, 9.0, 7.0, 9.0, 5.0, 10.0, 1.0, 9.0, 9.0, 9.0, 2.0, 1.0, 8.0, 10.0, 9.0, 1.0, 2.0, 9.0, 10.0, 8.0, 4.0, 3.0, 8.0, 8.0, 3.0, 6.0, 8.0, 5.0, 1.0, 5.0, 8.0, 6.0, 4.0, 8.0, 6.0, 4.0, 5.0, 5.0]
global b_y = 10
global p = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
global q = [0.827, 0.8, 0.977, 0.935, 0.996, 0.924, 0.943, 0.913, 0.908, 0.933, 0.907, 0.921, 0.886, 0.879, 0.841, 0.904, 0.917, 0.889, 0.936, 0.895, 0.878, 0.851, 0.981, 0.907, 0.831, 0.834, 0.964, 0.848, 0.891, 0.889, 0.838, 0.89, 0.86, 0.831, 0.898, 0.836, 0.934, 0.955, 0.824, 0.935, 0.985, 0.989, 0.857, 0.902, 0.97, 0.927, 0.85, 0.866, 0.868, 0.824, 0.992, 0.81, 0.946, 0.835, 0.951, 0.993, 0.919, 0.802, 0.851, 0.984, 0.982, 0.805, 0.952, 0.978, 0.92, 0.914, 0.82, 0.929, 0.984, 0.944, 0.997, 0.945, 0.891, 0.852, 0.992, 0.919, 0.982, 0.89, 0.875, 0.86, 0.929, 0.847, 0.971, 0.95, 0.87, 0.924, 0.979, 0.878, 0.911, 0.997, 0.867, 0.87, 0.871, 0.909, 0.927, 0.963, 0.973, 0.883, 0.93, 0.902, 0.962, 0.97, 0.906, 0.918, 0.938, 0.867, 0.885, 0.808, 0.95, 0.875, 0.919, 0.898, 0.829, 0.846, 0.821, 0.804, 0.942, 0.928, 0.928, 0.914, 0.92, 0.983, 0.957, 0.993, 0.931, 0.942, 0.802, 0.916, 0.935, 0.923, 0.821, 0.918, 0.952, 0.844, 0.936, 0.929, 0.872, 0.874, 0.991, 0.863, 0.994, 0.978, 0.917, 0.885, 0.883, 0.811, 0.809, 0.895, 0.869, 0.994, 0.902, 0.978, 0.998, 0.841, 0.891, 0.849, 0.865, 0.894]
global origin = 1
global destination = 40