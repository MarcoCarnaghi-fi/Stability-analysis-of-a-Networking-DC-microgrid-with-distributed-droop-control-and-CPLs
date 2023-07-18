% THIS SCRIPT CORRESPOND TO THE RESULTS A

% Script for the analysis of a 3-unit Networking DCMG without considering Interconnection impedances

clc 
clear
%------------------------------------------------------------
%   Parameter definition
%------------------------------------------------------------

% General
V_ref = 800;
n = 3;
tau = 0.003/5;

% Passive filter
    modC = 1.0;
    C1 = 0.33e-3    * modC ;
    C2 = 0.5e-3    * modC ;
    C3 = 1e-3      * modC ;
    
    modL = 1.0;
    L1 = 1.0e-3    * modL  ;
    L2 = 0.85e-3   * modL   ;
    L3 = 1.15e-3   * modL  ;
    
    rs1 = 0.45       ;
    rs2 = 0.525      ;
    rs3 = 0.5        ;
    
% Interconnection
    r_jk = 2          ;
    g12 = 0.5         ;
    g13 = 0           ;
    g14 = 0.5         ;
    g15 = 0           ;
    g16 = 0           ;
    g23 = 0.5         ;
    g24 = 0.5         ;
    g25 = 0           ;
    g26 = 0           ;
    g34 = 0           ;
    g35 = 0           ;
    g36 = 0.5         ;
    g45 = 0.5         ;
    g46 = 0           ;
    g56 = 0.5         ;
    
    aij = 1;
    a12 = aij;
    a13 = aij;
    a23 = aij;
    
    A_consensus = [0, a12, a13;...
                   0, 0  , a23;...
                   0  , 0  , 0 ];
   A_consensus = A_consensus + A_consensus';
   
   LA = diag(sum(A_consensus)) - A_consensus;
   
% Control

    k_iv_1  =   10        ; 
    k_iv_2  =   10        ;
    k_iv_3  =   10          ;
    
    k_pv_1  =   0.01             ; 
    k_pv_2  =   0.01           ;
    k_pv_3  =   0.01           ;
    
    k_ii_1  =   10              ;
    k_ii_2  =   10          ;
    k_ii_3  =   10          ;
     
    
    k_pi_1  =   0.01          ;
    k_pi_2  =   0.01          ;
    k_pi_3  =   0.01          ;
    
   
    
    I_rate  =   [10, 20, 30]    ;
    k_1     =   1/I_rate(1)     ;
    k_2     =   1/I_rate(2)     ;
    k_3     =   1/I_rate(3)     ;
    
    r_d0    =   0.3                 ;
    kd      =   [1, 2, 3]  ;
    rd1     =   r_d0/kd(1)          ;
    rd2     =   r_d0/kd(2)          ;
    rd3     =   r_d0/kd(3)          ;
    
  
%-----------------------------------------------------------------
% MATRIX POLYTOPIC GENERATION        
%-----------------------------------------------------------------
        vertex          = 2^n;
        P_CPL_max       = 160e3;
        P_CPL_min       = 10e3;
        P_DER_max       = 0e3;
        P_DER_min       = 0e3;
        %P_CPL_sup = P_CPL_max*[1 , 0.25  , 0.1 ]; %Unbalanced
        P_CPL_sup = P_CPL_max*[0.1 , 0.25  , 1 ]; %Balanced
        
        P_CPL_inf = P_CPL_min*[1 ,1 ,1];
        P_DER_sup = P_DER_max*[1 , 0   , 0.5  ];
        P_DER_inf = P_DER_min*[1 , 0   , 2    ];
        
        P_load_sup  = P_CPL_sup - P_DER_inf ;
        P_load_inf  = P_CPL_inf - P_DER_sup ;

         p_aux = [P_load_sup;P_load_inf];
         p_load_vertex=[];

         for i1=1:2
            for i2=1:2
                for i3=1:2
                     p_load_vertex = [p_load_vertex;...
                                      p_aux(i1,1) , p_aux(i2,2) , p_aux(i3,3)]; 
                
                        
                end
            end
         end
%-----------------------------------------------------------------
% SYTEM MATRIX GENERATION        
%-----------------------------------------------------------------    
    vertex = 2^n;
    A_glob = zeros(vertex,4*n-(n-1),4*n-(n-1)); 
    eigenvalues = [];
    
    
    for i=1:vertex
    
    
         
         G_CPL = -p_load_vertex(i,:)/(V_ref^2);
         G_CPL = diag(G_CPL);
         
         %Inductance
         IND    = [L1,L2,L3];
         IND    = IND(1:n);
         m1_L    = diag(IND.^(-1));
         
         %Capacitance
         CAP        = 1.0* [C1, C2, C3]  ;
         CAP        = CAP(1:n)                  ;
         m1_C       = diag(CAP.^(-1));
         m1_C2      = m1_C * ones(n,n);
        
         %Droop resistance
         RD0 = [rd1,rd2,rd3] * diag(ones(1,n));
         %DELTA_RD_EQ = zeros(n,n);
         
         %Output resistance
         R_SK    = [rs1,rs2,rs3];
         R_SK    = diag(R_SK(1:n));
         
         %Interconnection resistance
         G_JK    = [0   ,g12, g13;...
                    g12 ,0  , g23;...
                    g13 ,g23, 0  ];
         G_JK   = -G_JK(1:n,1:n);
         G_JK   = G_JK + diag(-1*sum(G_JK)); 
         mLE    = G_JK;
         
         %I linealization
         I_EQ 		= zeros(1,n)    			;%get from sim
        
         for k=1:n
            DELTA_RD_EQ(k)	=  (0.1) .* (kd(k)^(-1));         
            I_EQ(k) 		=  sum(p_load_vertex(i,:)) * (kd(k)/(sum(kd)));
            I_EQ(k)         =  I_EQ(k) ./ (V_ref);
         end
         I_EQ = diag(I_EQ);
         
         ALPHA = zeros(1,n);
            for k=1:n
                %ALPHA(k) = I_rate(k)/(sum(I_rate.*NEIGHBOURHOOD_MATRIX(k,:)));
                ALPHA(k) = I_rate(k)/(sum(I_rate));
            end
         ALPHA = diag(ALPHA);
         
         
         DELTA_RD_EQ = diag(DELTA_RD_EQ);
         
         %Control Parameter
         mK_II   = diag([k_ii_1,k_ii_2,k_ii_3]);
         mK_PI   = diag([k_pi_1,k_pi_2,k_pi_3]);
         mK_IV   = diag([k_iv_1,k_iv_2,k_iv_3]);
         mK_PV   = diag([k_pv_1,k_pv_2,k_pv_3]);
         
         %zero matrix
         cero = zeros(n,n);
         
         %A1 are the rows related to the i variable
         A1 = [-m1_L*(RD0+DELTA_RD_EQ+R_SK),    -m1_L ,  -m1_L*I_EQ     ,   m1_L];
         
         %A2 are the rows related to the v variable
         A2 = [m1_C , -m1_C2 .*( G_CPL) , cero, cero];
         
         %A3 are the rows related to  DeltaRd
         A3     = [mK_II,  mK_II* ALPHA * (ones(n,n) * G_CPL), cero , cero];
         A3_aux2= mK_PI * A1;
         A3 = A3 + A3_aux2;
         
         
         %A4 are the rows related to Deltav
         A4 = [cero, -1/n .* mK_IV *ones(n,n) , cero ,cero];
         A4_aux = -1/n.* mK_PV *(ones(n,1) * sum(A2)); 
         A4 = A4 + A4_aux;
                  
         
         A = [A1;A2;A3;A4];
         
         
         A(:,n+1) = A(:,n+1) + A(:,n+2)+ A(:,n+3); %All voltage are the same variable in a SinglePCC scenario
         A(n+1,:) = A(n+1,:) + A(n+2,:)+ A(n+3,:); 
         A = [A(:,1:n+1),A(:,n+4:end)];
         
         A = [A(1:n+1,:);A(n+4:end,:)];
         eigenvalues = [eigenvalues , eig(A)];

         A_glob(i,:,:) = A;
         
    end  
      

%-----------------------------------------------------------------
% STORE EIGENVALUES     
%-----------------------------------------------------------------  
eig_single = [eigenvalues(1:(end-3),1),eigenvalues(1:(end-3),end)]
csvwrite('eig_single.csv',eig_single);
%-----------------------------------------------------------------
% SIMULATION TIMES  
%-----------------------------------------------------------------  

    
    tau = 0.003/5 ;
    time1 = 0.5;
    time2 = 1.5;
    time3 = 5;