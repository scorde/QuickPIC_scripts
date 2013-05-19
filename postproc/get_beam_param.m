function beam_param = get_beam_param( SUB_BEAM, npart )


beam_param.Fraction = size(SUB_BEAM,1)/npart;
for i=1:size(SUB_BEAM,2)
    tmp(:,:) = SUB_BEAM(:,i,[1,4]);
    Cx = cov(tmp);
    tmp(:,:) = SUB_BEAM(:,i,[2,5]);
    Cy = cov(tmp);
    beam_param.eps_x(i) = sqrt(det(Cx));
    beam_param.eps_y(i) = sqrt(det(Cy));
    beam_param.sig_x(i) = sqrt(Cx(1,1));
    beam_param.sig_xp(i) = sqrt(Cx(2,2));
    beam_param.sig_y(i) = sqrt(Cy(1,1));
    beam_param.sig_yp(i) = sqrt(Cy(2,2));
    beam_param.beta_x(i) = Cx(1,1) / beam_param.eps_x(i);
    beam_param.alpha_x(i) = -Cx(1,2) / beam_param.eps_x(i);
    beam_param.beta_y(i) = Cy(1,1) / beam_param.eps_y(i);
    beam_param.alpha_y(i) = -Cy(1,2) / beam_param.eps_y(i); 
    beam_param.eps_nx(i) = mean(1e3*SUB_BEAM(:,i,6)/0.511) * beam_param.eps_x(i);
    beam_param.eps_ny(i) = mean(1e3*SUB_BEAM(:,i,6)/0.511) * beam_param.eps_y(i);
end

end

