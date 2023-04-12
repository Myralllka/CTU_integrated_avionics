function [Cbn] = norm_DCM(DCMorig)
    % function to orthogonalize and normalize DCM matrix
    % DCM = DCMorig;
    % DCM = DCM + 0.5*(eye(3) - DCM*DCM')*DCM;    % orthogonalization;
    % 
    % DCM(1,:)=DCM(1,:)/VectorNormC(DCM(1,:),3);  % normalization of DCM by rows
    % DCM(2,:)=DCM(2,:)/VectorNormC(DCM(2,:),3);
    % DCM(3,:)=DCM(3,:)/VectorNormC(DCM(3,:),3);

    Cbn = DCMorig;
    delta_12 = Cbn(1,:)*Cbn(2,:)';
    Cbn(1,:) = Cbn(1,:) - 1/2*delta_12*Cbn(2,:);
    Cbn(2,:) = Cbn(2,:) - 1/2*delta_12*Cbn(1,:);
    Cbn(3,:) = cross(Cbn(1,:),Cbn(2,:));
    Cbn(1,:) = Cbn(1,:)./norm(Cbn(1,:));
    Cbn(2,:) = Cbn(2,:)./norm(Cbn(2,:));
    Cbn(3,:) = Cbn(3,:)./norm(Cbn(3,:));

end