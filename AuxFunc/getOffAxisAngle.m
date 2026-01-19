function varPhi = getOffAxisAngle(groupAzEl,userAzEl)
    [~,nUser,nGroup] = size(userAzEl);
    % PhiTheta为两行矩阵，第一行代表方位角phi，第二行代表俯仰角theta，第i列代表用户i
    varPhi = zeros(nGroup,nUser,nGroup);
    for i = 1:nGroup
        for n = 1:nUser
            temp = repmat(userAzEl(:,n,i),1,nGroup);
            varPhi(:,n,i) = calculateTheta(groupAzEl,temp);
        end
    end
end

function theta = calculateTheta(PhiTheta1,PhiTheta2)
    theta = acosd(cosd(PhiTheta1(1,:)).*cosd(PhiTheta2(1,:))+sind(PhiTheta1(1,:)).*sind(PhiTheta2(1,:)).*cosd(PhiTheta1(2,:)-PhiTheta2(2,:)));
end