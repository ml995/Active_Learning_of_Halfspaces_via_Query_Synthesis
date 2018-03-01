function [alpha_1, alpha_2, tot_it] = noisy_find_coeff_o2(e_1, e_2, orackle, epsilon, noise, beta)
format long

an = [dot(orackle,e_1),dot(orackle,e_2)];
an = an/norm(an);

%temp = [0.5,0.5];
temp = rand(1,2);
temp = temp/norm(temp);
d = noisy_compute_sign(orackle, temp(1,2)*e_1-temp(1,1)*e_2,noise);

% l is a linked list. Each node in l represents an interval on the unit
% circle and has six entries. The first two entries denote the coordinate
% of the point that the interval begins with and the 3rd and 4th entries
% denote the coordinate of the point that the interval ends up with.
% The 5th entry denotes the probability density on the interval.
% The 6th entry denotes the probability mass of the interval.
l = java.util.LinkedList();
temp = d * temp;
l.addLast([-temp,temp,1-beta,2*(1-beta)]);
l.addLast([temp,-temp,beta,2*beta]);
n = 2;

for tot_it = 1:1100
    it = l.listIterator();
    s = 0;
    % Here we start with a node in the linked list. n is the total #
    % of the nodes in the linked list. We compute the total probability
    % mass of consecutive n/2 nodes/intervals. Then we substract the
    % probability mass of the first node and add the probability mass
    % of the next node, and so on, in order to find out when the total
    % probability mass of consecutive n/2 intervals change from >1/2 to <1/2
    % or from <1/2 to >1/2. Then we can find the interval where we will
    % find the zero point and that we will divide. Along the diameter that
    % passes the zero point, we cut the unit circle into two halves and
    % each half contains an exact probability mass of 1/2 and now there
    % are n+2 intervals.
    for ii = 1:(n/2)
        C = it.next();
        s = s + C(5,1);
    end

    sgn_prev = sign(s-1/2);
    it2 = l.listIterator();
    
    for ii = 1:(n/2)
        C2 = it2.next();
        C1 = it.next();
        old_s = s;
        % We substract the first node and add the next node.
        s = s - C2(5,1) + C1(5,1);
        sgn_now = sign(s-1/2);
        if (sgn_now ~= sgn_prev)
            break;
        end
        sgn_prev = sgn_now;
    end
    
    % We calculate where the zero point resides.
    x = (1/2 - old_s)/(C1(5,1) - C2(5,1));
    z = dot(C1(3:4,1),C1(1:2,1));
    if z > 1
        z = 1;
    end
    if z < -1
        z = -1;
    end
    angle = acos(z) * x;
    bisec = [cos(angle),-sin(angle);sin(angle),cos(angle)]*C1(1:2,1);
    d = noisy_compute_sign(orackle, bisec(2,1)*e_1-bisec(1,1)*e_2,noise);
    
    if d > 0
        left = 2*(1-beta); right = 2*beta;
    else
        left = 2*beta; right = 2*(1-beta);
    end
    
    C = C2;
    C2 = [C2(1:2,1);-bisec;x*C(5,1)*right;C(6,1)*right];
    it2.set(C2);
    C_temp = [-bisec;C(3:4,1);(1-x)*C(5,1)*left;C(6,1)*left];
    it2.add(C_temp);
    
    for ii = 1:(n-2)/2
        C = it2.next;
        it2.set([C(1:4,1); C(5,1) * left;C(6,1)*left]);
    end
    it2.next();
    C = C1;
    C1 = [C1(1:2,1);bisec;x*C(5,1)*left;C(6,1)*left];
    it2.set(C1);
    C_temp = [bisec;C(3:4,1);(1-x)*C(5,1)*right;C(6,1)*right];
    it2.add(C_temp);
    
    for ii = 1:(n-2)/2
        if it2.hasNext == 0
            it2 = l.listIterator;
        end
        C = it2.next;
        it2.set([C(1:4,1); C(5,1) * right;C(6,1)*right]);
    end
    n = l.size();
    it = l.listIterator();
    C_max = [0 0 0 0 0 0]';
    % Here we find the interal with highest probability mass.
    while it.hasNext()
        C_temp = it.next;
        if C_temp(5,1)>C_max(5,1)%highest probability density
            C_max = C_temp;
            
        end
    end
    
    temp = (C_max(1:2,1)+C_max(3:4,1))/2;
    temp = temp/norm(temp);
    
    if norm(temp'-an)<epsilon
        break;
    end
end

alpha_1 = temp(1,1);
alpha_2 = temp(2,1);





