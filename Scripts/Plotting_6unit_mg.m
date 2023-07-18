% Script to plot time curves of Voltage and Current from Simulink simulation


    time_sa = V_MG_sa.Time;
    figure(1)
    plot(time_sa,   v_b1_sa.Data,'Linewidth',3.0);
    hold on
    plot(time_sa,   v_b2_sa.Data, 'Linewidth',3.0);
    plot(time_sa,   v_b3_sa.Data, 'Linewidth',3.0);
    plot(time_sa,   v_b4_sa.Data, 'Linewidth',3.0);
    plot(time_sa,   v_b5_sa.Data, 'Linewidth',3.0);
    plot(time_sa,   v_b6_sa.Data, 'Linewidth',3.0);
    plot(time_sa,   V_MG_sa.Data, 'Linewidth',3.0);
    %plot(time_sa,V_MG_sa.Data, 'k','Linewidth',2.0);
    grid


    title('Voltages')

    voltages_6mg(:,1) = time_sa;
    voltages_6mg(:,2) = v_b1_sa.Data;
    voltages_6mg(:,3) = v_b2_sa.Data;
    voltages_6mg(:,4) = v_b3_sa.Data;
    voltages_6mg(:,5) = v_b4_sa.Data;
    voltages_6mg(:,6) = v_b5_sa.Data;
    voltages_6mg(:,7) = v_b6_sa.Data;
    voltages_6mg(:,8) = V_MG_sa.Data;

    csvwrite('voltages_6mg.csv',voltages_6mg);
    

    figure(2)
    plot(time_sa,   i_1_sa.Data, 'b','Linewidth',3.0);
    hold on
    plot(time_sa,   i_2_sa.Data, 'r','Linewidth',3.0);
    plot(time_sa,   i_3_sa.Data, 'k','Linewidth',3.0);
    plot(time_sa,   i_4_sa.Data, 'b','Linewidth',3.0);
    plot(time_sa,   i_5_sa.Data, 'r','Linewidth',3.0);
    plot(time_sa,   i_6_sa.Data, 'k','Linewidth',3.0);
    grid
    title('Current')

    currents_6mg(:,1) = time_sa;
    currents_6mg(:,2) = i_1_sa.Data;
    currents_6mg(:,3) = i_2_sa.Data;
    currents_6mg(:,4) = i_3_sa.Data;
    currents_6mg(:,5) = i_4_sa.Data;
    currents_6mg(:,6) = i_5_sa.Data;
    currents_6mg(:,7) = i_6_sa.Data;


    csvwrite('currents_6mg.csv',currents_6mg);