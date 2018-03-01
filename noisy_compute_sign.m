function x = noisy_compute_sign (y, orackle, noise)
    x = sign(dot(y,orackle));
    if rand() < noise
        x = - x;
    end
end