function [features, distributions]=compute_features(t,X,pts)
% given a list of (t,x) pts that indicate periods in the traces X, compute
% per period feature distributions and their summary statistics 