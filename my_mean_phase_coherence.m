

function mpc = my_mean_phase_coherence(data1,data2)
        phase1=angle(hilbert(data1));
        phase2=angle(hilbert(data2));
        ph_diff=phase2-phase1;
        mpc = sqrt(mean(cos(ph_diff)).^2 + mean(sin(ph_diff)).^2);
end

