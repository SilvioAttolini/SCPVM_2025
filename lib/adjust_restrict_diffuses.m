function [c_scaled, d_scaled] = adjust_restrict_diffuses(a, b, c, d)

    power_a = mean(a.^2);
    power_b = mean(b.^2);
    power_c = mean(c.^2);
    power_d = mean(d.^2);

    target_power = (power_a + power_b) / 2;

    % Scale c and d to match the target power
    scaling_factor_c = sqrt(target_power / power_c);
    scaling_factor_d = sqrt(target_power / power_d);
    
    % Apply the scaling factors
    c_scaled = c * scaling_factor_c;
    d_scaled = d * scaling_factor_d;

end
