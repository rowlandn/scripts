function [m9_init_idx, d_plot] = findDip(s)

[max_val, max_idx] = calcMax(s)
[min_val, min_idx] = calcMinVm(s, max_idx)

[m9_init_idx, d_plot] = calcInitVmV3hKpTinterp(s, max_idx, min_idx, 5, 50, 1);