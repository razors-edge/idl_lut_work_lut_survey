pro angular_distance_array,alpha1_d,beta1_d,alpha2_d,beta2_d,d_deg

    ind = 3.1415926 / 180.0
    alpha1_radian=alpha1_d  * ind
    beta1_radian=beta1_d * ind
    alpha2_radian=alpha2_d * ind 
    beta2_radian=beta2_d * ind
    
    delta_alpha = alpha1_radian - alpha2_radian
    delta_beta = beta1_radian - beta2_radian
    hav_delta_alpha = sin(delta_alpha/2.)*sin(delta_alpha/2.) 
    hav_delta_beta = sin(delta_beta/2.)*sin(delta_beta/2.) 
    hav_d = hav_delta_beta + (cos(beta1_radian)*cos(beta2_radian)*hav_delta_alpha)
    d_radian = acos(1 - (hav_d *2.))
    d_deg = d_radian / ind
    
    end
