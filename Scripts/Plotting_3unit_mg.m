% Script to plot time curves of Voltage and Current from Simulink simulation


%% 
%-----------------------------------------------------------------
% UNBALANCED 
%-----------------------------------------------------------------  
    time_sa = V_MG_sa.Time;
    figure(1)
    plot(time_sa,v_b1_sa.Data,'Linewidth',3.0);
    hold on
    plot(time_sa,v_b2_sa.Data, 'Linewidth',3.0);
    plot(time_sa, v_b3_sa.Data, 'Linewidth',3.0);
    plot(time_sa,V_MG_sa.Data, 'Linewidth',3.0);
    %plot(time_sa,V_MG_sa.Data, 'k','Linewidth',2.0);
    grid


    title('Voltages')

    voltages_unbalanced(:,1) = time_sa;
    voltages_unbalanced(:,2) = v_b1_sa.Data;
    voltages_unbalanced(:,3) = v_b2_sa.Data;
    voltages_unbalanced(:,4) = v_b3_sa.Data;
    voltages_unbalanced(:,5) = V_MG_sa.Data;

    csvwrite('voltages_unbalanced.csv',voltages_unbalanced);
    

    figure(2)
    plot(time_sa,i_1_sa.Data, 'b','Linewidth',3.0);
    hold on
    plot(time_sa,i_2_sa.Data, 'r','Linewidth',3.0);
    plot(time_sa, i_3_sa.Data, 'k','Linewidth',3.0);
    grid
    title('Current')

    currents_unbalanced(:,1) = time_sa;
    currents_unbalanced(:,2) = i_1_sa.Data;
    currents_unbalanced(:,3) = i_2_sa.Data;
    currents_unbalanced(:,4) = i_3_sa.Data;

    csvwrite('currents_unbalanced.csv',currents_unbalanced);

%% 
%-----------------------------------------------------------------
% BALANCED 
%-----------------------------------------------------------------    

    time_sa = V_MG_sa.Time;
    figure(1)
    plot(time_sa,v_b1_sa.Data,'Linewidth',3.0);
    hold on
    plot(time_sa,v_b2_sa.Data, 'Linewidth',3.0);
    plot(time_sa, v_b3_sa.Data, 'Linewidth',3.0);
    plot(time_sa,V_MG_sa.Data, 'Linewidth',3.0);
    %plot(time_sa,V_MG_sa.Data, 'k','Linewidth',2.0);
    grid


    title('Voltages')

    voltages_balanced(:,1) = time_sa;
    voltages_balanced(:,2) = v_b1_sa.Data;
    voltages_balanced(:,3) = v_b2_sa.Data;
    voltages_balanced(:,4) = v_b3_sa.Data;
    voltages_balanced(:,5) = V_MG_sa.Data;

    csvwrite('voltages_balanced.csv',voltages_balanced);
    

    figure(2)
    plot(time_sa,i_1_sa.Data, 'b','Linewidth',3.0);
    hold on
    plot(time_sa,i_2_sa.Data, 'r','Linewidth',3.0);
    plot(time_sa, i_3_sa.Data, 'k','Linewidth',3.0);
    grid
    title('Current')

    currents_balanced(:,1) = time_sa;
    currents_balanced(:,2) = i_1_sa.Data;
    currents_balanced(:,3) = i_2_sa.Data;
    currents_balanced(:,4) = i_3_sa.Data;

    csvwrite('currents_balanced.csv',currents_balanced);

%% 
%-----------------------------------------------------------------
% SINGLEPCC 
%-----------------------------------------------------------------  

    time_sa = V_MG_sa.Time;
    figure(1)
    plot(time_sa,V_MG_sa.Data, 'Linewidth',3.0);
    grid
    title('Voltages')

    voltages_single(:,1) = time_sa;
    voltages_single(:,2) = V_MG_sa.Data;

    csvwrite('voltages_single.csv',voltages_single);
    

    figure(2)
    plot(time_sa,i_1_sa.Data, 'b','Linewidth',3.0);
    hold on
    plot(time_sa,i_2_sa.Data, 'r','Linewidth',3.0);
    plot(time_sa, i_3_sa.Data, 'k','Linewidth',3.0);
    grid
    title('Current')

    currents_single(:,1) = time_sa;
    currents_single(:,2) = i_1_sa.Data;
    currents_single(:,3) = i_2_sa.Data;
    currents_single(:,4) = i_3_sa.Data;

    csvwrite('currents_single.csv',currents_single);
    
    
    %% EIGENVALUES 3-UNIT DCMG UNBALANCED

eig_unb = [eigenvalues(:,1),eigenvalues(:,end)]