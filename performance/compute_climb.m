function climb = compute_climb(ExcessP_W, W, V)
% Computes rate of climb and climb angle

climb.ROC_mps = ExcessP_W / W;
ROC_over_V = climb.ROC_mps / V;
ROC_over_V(ROC_over_V > 1) = 1;
ROC_over_V(ROC_over_V < -1) = -1;
climb.gamma_deg = rad2deg(asin(ROC_over_V));

end
