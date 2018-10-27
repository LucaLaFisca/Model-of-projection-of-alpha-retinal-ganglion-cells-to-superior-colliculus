%% From Stimulus to SC Spiking Rates

%% ALPHA RGCs
%%
%% struct to save all results
ALL_RGC_RFmap = struct('session',1, 'Rates',zeros(500));
IDX_RGC_RFmap = 0;
ALL_STIM = struct('session',1, 'trial',1, 'stim',zeros(500));
IDX_STIM = 0;
ALL_RGC_RETSTIM = struct('session',1, 'Rates',zeros(500));
IDX_RGC_RETSTIM = 0;
ALL_WEIGHTS = struct('session',1, 'channel',1, 'gain',1, 't_shift',40, 'spont',40, 'weights',zeros(7));
IDX_WEIGHTS = 0;
ALL_RFmap_RATES = struct('session',1, 'channel',1, 'rates',zeros(500));
IDX_RFmap = 0;
ALL_wrong_WEIGHTS = struct('session',1, 'channel',1, 'gain',1, 't_shift',40, 'spont',40, 'weights',zeros(7));
IDX_wrong_WEIGHTS = 0;
ALL_wrong_RFmap_RATES = struct('session',1, 'channel',1, 'rates',zeros(500));
IDX_wrong_RFmap = 0;
ALL_RETSTIM_RATES = struct('session',1, 'channel',1, 'rates',zeros(500));
IDX_RETSTIM = 0;

%%
for session = 1:5   %5 different locations recorded
    disp('session')
    disp(session)
    close all
    
    %% SPATIO-TEMPORAL PARAMETERS
    dx = .7;    %degree
    frame_rate = 60*3;   %Hz
    dt = 1/frame_rate * 1000; %ms
    tau_vect = 0:dt:500; %ms %Span of temporal kernel
    
    %% Run settings Ephys
    global Set
    run runSettingsEphys

    screenx = Set.Screenx * Set.DegPerPix;  %screen width in degree
    screeny = Set.Screeny * Set.DegPerPix;  %screen height in degree
    mouseposx = (Set.ScreenWidth-Set.mouseposx) * (Set.Screenx/Set.ScreenWidth);  %from cm to px
    mouseposx = mouseposx * Set.DegPerPix; %from px to degree
    mouseposy = Set.mouseposy * (Set.Screeny/Set.ScreenHeight);
    mouseposy = mouseposy * Set.DegPerPix;

    %%Load the RetStim log
    % load('Test_20180912_B3.mat');
    clear RetstimMat
    my_path = ['Recordings\Session-', int2str(session), '\RetStimMat.mat'];
    load(my_path);
    My_stim = RetstimMat;
    My_stim.centerx = My_stim.centerx * Set.DegPerPix; %converts to degree
    My_stim.centery = My_stim.centery * Set.DegPerPix;
    stim_centerx = My_stim.centerx(1); %position of RF from the center of the screen
    stim_centery = My_stim.centery(1);
    %mouse position considering the center of RF as the origin
    mouseposx = -stim_centerx - (screenx/2 - mouseposx);
    mouseposy = -stim_centery - (screeny/2 - mouseposy);

    My_stim.centerx = (My_stim.centerx - stim_centerx); %position from center of RF
    My_stim.centery = (My_stim.centery - stim_centery);

    My_stim.size = My_stim.size * Set.DegPerPix;
    RF = max(My_stim.size)/3*2;

    My_stim.sprite(My_stim.sprite == 2) = -1;   %sprite 2 = black dot (cf.RetStim)

    My_stim.time = My_stim.time * 1000; %Converts to ms

    %%Load RetStim MUA (spike rates)
    my_path = ['Recordings\Session-', int2str(session), '\RetStimMUA.mat'];
    RetstimMUA = cell2mat(struct2cell(load(my_path)));
%     smooth_mua1 = smooth(RetstimMUA(:),30);
%     smooth_mua2 = smooth(smooth_mua1,30);
%     RetstimMUA = reshape(smooth_mua2, size(RetstimMUA,1),size(RetstimMUA,2),size(RetstimMUA,3));

    %%Load Retstim px (time vector)
    my_path = ['Recordings\Session-', int2str(session), '\RetStimPx.mat'];
    RetstimPx = cell2mat(struct2cell(load(my_path)));

    %%Load RFmap log
    % my_path = input('Enter path of "RFmapMat" data (from "recordings" directory)\n' ,'s');
    % Mat = cell2mat(struct2cell(load(my_path)));
    my_path = ['Recordings\Session-', int2str(session), '\RFmapMat.mat'];
    RFmapMat = cell2mat(struct2cell(load(my_path)));

    %%Load Multi-chanel spike rates
    % my_path = input('Enter path of "RFmapMUA" data (from "recordings" directory)\n' ,'s');
    % Mua = cell2mat(struct2cell(load(my_path)));
    my_path = ['Recordings\Session-', int2str(session), '\RFmapMUA.mat'];
    RFmapMUA = cell2mat(struct2cell(load(my_path)));
%     test = cell2mat(struct2cell(load(my_path)));
%     my_smooth = smooth(test(:),120);
%     MY_SMOOTH = reshape(my_smooth, size(test,1),size(test,2),size(test,3));

%     smooth_mua1 = smooth(RFmapMUA(:),30);
% %     SMOOTH_MUA1 = reshape(smooth_mua1, size(RFmapMUA,1),size(RFmapMUA,2),size(RFmapMUA,3));
%     smooth_mua2 = smooth(smooth_mua1,30);
%     RFmapMUA = reshape(smooth_mua2, size(RFmapMUA,1),size(RFmapMUA,2),size(RFmapMUA,3));

%     smooth_mua3 = smooth(RFmapMUA(:),60);
%     SMOOTH_MUA3 = reshape(smooth_mua3, size(RFmapMUA,1),size(RFmapMUA,2),size(RFmapMUA,3));    
%     figure; hold on; plot(px, RFmapMUA(:,20,20),'green'); ...
%         plot(px, SMOOTH_MUA1(:,20,20),'red');plot(px, SMOOTH_MUA2(:,20,20),'blue');...
%         plot(px, SMOOTH_MUA3(:,20,20),'black')

    %%Load time vector
    % my_path = input('Enter path of "RFmapPx" data (from "recordings" directory)\n' ,'s');
    % px = cell2mat(struct2cell(load(my_path)));
    my_path = ['Recordings\Session-', int2str(session), '\RFmapPx.mat'];
    px = cell2mat(struct2cell(load(my_path)));
    RFmapPx = cell2mat(struct2cell(load(my_path)));
    
    %% Define the plane (respecting to the screen size)
    % reference = getGlobalDataSpatial;
    % x_spatial = reference(1,:);
    % y_spatial = reference(2,:);
    % [~,idx_center] = max(y_spatial);
    % d_alpha = x_spatial(idx_center);
    % d_alpha = 4.5; %degrees | according to "smallest diam" (cf.RetStim)
    % d_alpha = 4.7255; %degrees | size to have 49 RGCs inside SC rf
    % alpha_size = sqrt(d_alpha^2/2);
    % nb_RGCs = ceil((RF/2)/alpha_size);
    nb_RGCs = 3; %nb of RGCs on each part from the center
    alpha_size = RF/(2*nb_RGCs+1)/2; %degrees | size to have 49 RGCs inside SC rf
    d_alpha = sqrt(2*alpha_size^2);


    nb_check = 12; %cf.RFmap
    max_x = max(max(RFmapMat(:,2:5:5*nb_check)));
    min_x = min(min(RFmapMat(:,2:5:5*nb_check)));
    max_y = max(max(RFmapMat(:,3:5:5*nb_check)));
    min_y = min(min(RFmapMat(:,3:5:5*nb_check)));
    limit_xy = 15; %degree | useless to compute further regarding the model (sigma_s ~ 5)
    x_vect = mouseposx+min_x-(nb_RGCs+0.5)*alpha_size:dx:mouseposx+max_x+(nb_RGCs+0.5)*alpha_size;
    x_vect(x_vect<-limit_xy-(nb_RGCs+0.5)*alpha_size | x_vect>limit_xy+(nb_RGCs+0.5)*alpha_size) = [];
    y_vect = mouseposy+min_y-(nb_RGCs+0.5)*alpha_size:dx:mouseposy+max_y+(nb_RGCs+0.5)*alpha_size;
    y_vect(y_vect<-limit_xy-(nb_RGCs+0.5)*alpha_size | y_vect>limit_xy+(nb_RGCs+0.5)*alpha_size) = [];

    [xx, yy] = meshgrid(x_vect,y_vect);
    yy = flipud(yy);    %To keep the same plane of coordinates


    %% ON-CENTER TRANSIENT
    disp('ON transient')
    kernel_sign = 1;    %ON cell

    %%RGCs PARAMETERS
    %rates from litterature (cf."2_Alpha RGCs types")
    spont_rate_ON_t = 2.5;  %Hz
    fig_x = 'ON_trans_spatial.fig';
    fig_t = 'ON_trans_temporal.fig';

    %%FIT SPIKING RATES CURVES (fig_x & fig_t)
    % x0 = [0.998383184969581 0.999989967185898 0.0366604841402724 ...
    %     0.0223233592203604 0.0334478163810849 0.0207797425569440 ...
    %     6774.53849218199 6704.28724900915 0.703647134882378 5.08956641055278];
    %     %dx=0.1 frame_rate=300
% 
%     x0 = [0.960030765820577,0.993403917721996,0.0463412536933561,...
%         0.0330567691558252,0.0266709898031438,0.0169944916254541,...
%         1199.23331298142,699.773775756765,0.616379830248048,5.77385043755806];
%         %previous dx=0.7 frame_rate=60
    
%    x0 = [0.960392364841583,0.999919268312941,0.0460176700618532,...
%        0.0305559074245359,0.0265214822808016,0.0143224081124178,...
%        1199.40614908337,699.714430776737,0.616438990102239,5.95420506801665];
%         %dx=0.7 frame_rate=60
        
%     x0 = [0.982940223113363,0.999975068358858,0.0473677484726349,...
%         0.0312041565128494,0.0278460392289396,0.0149620273275837,...
%         1199.78455022815,699.068672764979,0.618022277868822,5.62622499867810];
%         %dx=0.7 frame_rate=120

    x0 = [0.983833759740827,0.957869041533103,0.0477479765894229,...
        0.0301209235385489,0.0270354506607895,0.0138503374477721,...
        1170.47058205900,694.027540727645,0.783908060733699,4.73038834242088];
        %dx=0.7 frame_rate=180

    [B_cen, B_sur, alpha_cen, alpha_sur, beta_cen, beta_sur, ...
        A_c, A_s,sigma_c, sigma_s] = ...
        fitSpikingRates(dx,dt,x0,spont_rate_ON_t,fig_x,fig_t);

    x0 = [B_cen, B_sur, alpha_cen, alpha_sur, beta_cen, beta_sur, ...
        A_c, A_s,sigma_c, sigma_s];

    save('Modelization\x0_ON_t', 'x0')

    % dx = x_vect(2) - x_vect(1); %degree
    % dt = tau_vect(2) - tau_vect(1); %ms


    %%DESIGN THE SPATIAL FILTER (kernel D_x)
    % x_array = repmat(x_vect,size(x_vect,2),1);
    % D_x_cen = A_c*exp(-(x_array.^2 + (x_array').^2)/...
    %     (2*(sigma_c^2)))/(2*pi*sigma_c^2);
    % D_x_sur = A_s*exp(-(x_array.^2 + (x_array').^2)/...
    %     (2*(sigma_s^2)))/(2*pi*sigma_s^2);
    D_x_cen_ON_t = A_c*exp(-(xx.^2 + yy.^2)/...
        (2*(sigma_c^2)))/(2*pi*sigma_c^2);
    D_x_sur_ON_t = A_s*exp(-(xx.^2 + yy.^2)/...
        (2*(sigma_s^2)))/(2*pi*sigma_s^2);
    D_x = kernel_sign * (D_x_cen_ON_t - D_x_sur_ON_t);
    figure
    % [xx,yy] = meshgrid(x_vect, x_vect);
    surf(xx,yy,D_x)
    title('Kernel ON-trans')
    xlabel('x(degrees)')
    ylabel('y (degrees)')
    zlabel('D_x')

    %%DESIGN TEMPORAL FILTER
    D_t_cen_ON_t = alpha_cen^2*tau_vect.*exp(-alpha_cen*tau_vect) - ...
       B_cen*beta_cen^2*tau_vect.*exp(-beta_cen*tau_vect);
    D_t_sur_ON_t = alpha_sur^2*tau_vect.*exp(-alpha_sur*tau_vect) - ...
       B_sur*beta_sur^2*tau_vect.*exp(-beta_sur*tau_vect);
    figure
    subplot(2,1,1)
    plot(tau_vect,D_t_cen_ON_t)
    title('D-t-cen-ON-trans')
    xlabel('time(ms)')
    ylabel('D_t')
    subplot(2,1,2)
    plot(tau_vect,D_t_sur_ON_t)
    title('D-t-sur-ON-trans')
    xlabel('time(ms)')
    ylabel('D_t')

    %%DESIGN THE SPATIO-TEMPORAL FILTER (kernel D_x_t)
    D_x_t_ON_trans = zeros(length(y_vect),length(x_vect),length(tau_vect));
    for k = 1:length(tau_vect)
        D_x_t_ON_trans(:,:,k) = kernel_sign * (D_t_cen_ON_t(k)*D_x_cen_ON_t -...
            D_t_sur_ON_t(k)*D_x_sur_ON_t);
    end
    [txt,xxt] = meshgrid(tau_vect,x_vect);
    plot_D_x_t = permute(D_x_t_ON_trans, [2 3 1]);
    figure
    % surf(xx,yy,plot_D_x_t(:,:,round(length(x_vect)/2)))
    [~,idx_y] = (min(abs(yy(:,1))));
    surf(xxt,txt,plot_D_x_t(:,:,idx_y))
    title('D-x-t-ON-trans')
    xlabel('x(degrees)')
    ylabel('time(ms)')
    zlabel('D_x_t')

    D_x_t_ON_trans = reshape(D_x_t_ON_trans, 1,length(y_vect)*length(x_vect)*length(tau_vect));


    %%ON-CENTER SUSTAINABLE
    disp('ON sustainable')
    kernel_sign = 1;    %ON

    %%RGCs PARAMETERS
    %rates from litterature (cf."2_Alpha RGCs types")
    spont_rate_ON_s = 25.5;   %Hz
    fig_x = 'ON_sust_spatial.fig';
    fig_t = 'ON_sust_temporal.fig';

    %%FIT SPIKING RATES CURVES (fig_x & fig_t)
    % x0 = [0.900500363205036 0.998813135209984 0.0516207460405460 ...
    %     0.0348737296920597 0.0252571829875721 0.0190306432706683 ...
    %     636.171321972219 476.841664259293 0.789091277373766 3.82829119811768];
    %     %10 parameters : dx=0.1 frame_rate=300

%     x0 = [0.865265278186953,0.981342748838577,0.0479500480170850,...
%         0.0304942721522636,0.0240729178711564,0.0177958196462081,...
%         643.426965113774,476.032437542049,0.589653899853549,4.58541383990350];
%         %previous dx=0.7 frame_rate=60

%     x0 = [0.862449788544349,0.957402853892592,0.0483291875186938,...
%         0.0315214425839663,0.0242323172465637,0.0188187704327758,...
%         643.547041498049,476.001296450496,0.628899605852294,4.62314098935297];
%         %dx=0.7 frame_rate=60

%     x0 = [0.905278225894168,0.921163149182086,0.0471326834589211,...
%         0.0381646324689537,0.0273952740277818,0.0245423384390752,...
%         764.374343462780,472.695134982016,0.635770927384516,4.62654836694358];
%         %dx=0.7 frame_rate=120

    x0 = [0.903128477758515,0.876365544433131,0.0521736379567930,...
        0.0621141678870341,0.0303939361655133,0.0412778480986692,...
        785.809694386541,487.500331131079,0.723710236527917,3.75540234769865];
        %dx=0.7 frame_rate=180
        
    [B_cen, B_sur, alpha_cen, alpha_sur, beta_cen, beta_sur, ...
        A_c, A_s,sigma_c, sigma_s] = ...
        fitSpikingRates(dx,dt,x0,spont_rate_ON_s,fig_x,fig_t);

    x0 = [B_cen, B_sur, alpha_cen, alpha_sur, beta_cen, beta_sur, ...
        A_c, A_s,sigma_c, sigma_s];

    save('Modelization\x0_ON_s', 'x0')

    %%DESIGN THE SPATIAL FILTER (kernel D_x)
    D_x_cen_ON_s = A_c*exp(-(xx.^2 + yy.^2)/...
        (2*(sigma_c^2)))/(2*pi*sigma_c^2);
    D_x_sur_ON_s = A_s*exp(-(xx.^2 + yy.^2)/...
        (2*(sigma_s^2)))/(2*pi*sigma_s^2);
    D_x = kernel_sign * (D_x_cen_ON_s - D_x_sur_ON_s);
    figure
    surf(xx,yy,D_x)
    title('Kernel ON-sust')
    xlabel('x(degrees)')
    ylabel('y (degrees)')
    zlabel('D_x')

    %%DESIGN TEMPORAL FILTER
    D_t_cen_ON_s = alpha_cen^2*tau_vect.*exp(-alpha_cen*tau_vect) - ...
       B_cen*beta_cen^2*tau_vect.*exp(-beta_cen*tau_vect);
    D_t_sur_ON_s = alpha_sur^2*tau_vect.*exp(-alpha_sur*tau_vect) - ...
       B_sur*beta_sur^2*tau_vect.*exp(-beta_sur*tau_vect);
    figure
    subplot(2,1,1)
    plot(tau_vect,D_t_cen_ON_s)
    title('D-t-cen-ON-sust')
    xlabel('time(ms)')
    ylabel('D_t')
    subplot(2,1,2)
    plot(tau_vect,D_t_sur_ON_s)
    title('D-t-sur-ON-sust')
    xlabel('time(ms)')
    ylabel('D_t')

    %%DESIGN THE SPATIO-TEMPORAL FILTER (kernel D_x_t)
    D_x_t_ON_sust = zeros(length(y_vect),length(x_vect),length(tau_vect));
    for k = 1:length(tau_vect)
        D_x_t_ON_sust(:,:,k) = kernel_sign * (D_t_cen_ON_s(k)*D_x_cen_ON_s -...
            D_t_sur_ON_s(k)*D_x_sur_ON_s);
    end
    [txt,xxt] = meshgrid(tau_vect,x_vect);
    plot_D_x_t = permute(D_x_t_ON_sust, [2 3 1]);
    figure
    [~,idx_y] = (min(abs(yy(:,1))));
    surf(xxt,txt,plot_D_x_t(:,:,idx_y))
    title('D-x-t-ON-sust')
    xlabel('x(degrees)')
    ylabel('time(ms)')
    zlabel('D_x_t')

    D_x_t_ON_sust = reshape(D_x_t_ON_sust, 1,length(y_vect)*length(x_vect)*length(tau_vect));

    %%OFF-CENTER TRANSIENT
    disp('OFF transient')
    kernel_sign = -1;   %OFF

    %%RGCs PARAMETERS
    %rates from litterature (cf."2_Alpha RGCs types")
    spont_rate_OFF_t = 10;  %Hz
    fig_x = 'OFF_trans_spatial.fig';
    fig_t = 'OFF_trans_temporal.fig';

    %%FIT SPIKING RATES CURVES (fig_x & fig_t)
    % x0 = [0.978686899348437 0.977098675977756 0.0578456657116329 ...
    %     0.0386300187502710 0.0207298418461319 0.0151202395002490 ...
    %     530.373243904072 179.538219382925 0.555872095479203 3.38867060322331];
    %     %dx=0.1 frame_rate=300

%     x0 = [0.934140512628427,0.933214617006371,0.0519137164500029,...
%         0.0474466861988930,0.0209962114479477,0.0221716487925633,...
%         618.945861540076,187.262878792505,0.560382737835537,3.90293295139530];
%         %previous dx=0.7 frame_rate=60

%     x0 = [0.936959397195609,0.938200936711232,0.0507651017019744,...
%         0.0428879968774772,0.0212821362690652,0.0196808469510948,...
%         638.943016190080,184.415216733833,0.560299244200916,3.90450115082461];
%         %dx=0.7 frame_rate=60
        
%     x0 = [0.974857373827497,0.974450225306137,0.0495847441663187,...
%         0.0669881420152485,0.0241114080639234,0.0306797056450317,...
%         726.675660607332,210.726236993740,0.563318282845221,3.88185394118526];
%         %dx=0.7 frame_rate=120

    x0 = [0.980851254610049,0.948893856087863,0.0492902107067040,...
        0.0890399065634060,0.0244409146466022,0.0355175070148930,...
        726.280386285538,210.561781351905,0.526300775390483,4.47981799723477];
        %dx=0.7 frame_rate=180
        
    [B_cen, B_sur, alpha_cen, alpha_sur, beta_cen, beta_sur, ...
        A_c, A_s,sigma_c, sigma_s] = ...
        fitSpikingRates(dx,dt,x0,spont_rate_OFF_t,fig_x,fig_t);

    x0 = [B_cen, B_sur, alpha_cen, alpha_sur, beta_cen, beta_sur, ...
        A_c, A_s,sigma_c, sigma_s];

    save('Modelization\x0_OFF_t', 'x0')

    %%DESIGN THE SPATIAL FILTER (kernel D_x)
    D_x_cen_OFF_t = A_c*exp(-(xx.^2 + yy.^2)/...
        (2*(sigma_c^2)))/(2*pi*sigma_c^2);
    D_x_sur_OFF_t = A_s*exp(-(xx.^2 + yy.^2)/...
        (2*(sigma_s^2)))/(2*pi*sigma_s^2);
    D_x = kernel_sign * (D_x_cen_OFF_t - D_x_sur_OFF_t);
    figure
    surf(xx,yy,D_x)
    title('Kernel OFF-trans')
    xlabel('x(degrees)')
    ylabel('y (degrees)')
    zlabel('D_x')

    %%DESIGN TEMPORAL FILTER
    D_t_cen_OFF_t = alpha_cen^2*tau_vect.*exp(-alpha_cen*tau_vect) - ...
       B_cen*beta_cen^2*tau_vect.*exp(-beta_cen*tau_vect);
    D_t_sur_OFF_t = alpha_sur^2*tau_vect.*exp(-alpha_sur*tau_vect) - ...
       B_sur*beta_sur^2*tau_vect.*exp(-beta_sur*tau_vect);
    figure
    subplot(2,1,1)
    plot(tau_vect,D_t_cen_OFF_t)
    title('D-t-cen-OFF-trans')
    xlabel('time(ms)')
    ylabel('D_t')
    subplot(2,1,2)
    plot(tau_vect,D_t_sur_OFF_t)
    title('D-t-sur-OFF-trans')
    xlabel('time(ms)')
    ylabel('D_t')

    %%DESIGN THE SPATIO-TEMPORAL FILTER (kernel D_x_t)
    D_x_t_OFF_trans = zeros(length(y_vect),length(x_vect),length(tau_vect));
    for k = 1:length(tau_vect)
        D_x_t_OFF_trans(:,:,k) = kernel_sign * (D_t_cen_OFF_t(k)*D_x_cen_OFF_t -...
            D_t_sur_OFF_t(k)*D_x_sur_OFF_t);
    end
    [txt,xxt] = meshgrid(tau_vect,x_vect);
    plot_D_x_t = permute(D_x_t_OFF_trans, [2 3 1]);
    figure
    [~,idx_y] = (min(abs(yy(:,1))));
    surf(xxt,txt,plot_D_x_t(:,:,idx_y))
    title('D-x-t-OFF-trans')
    xlabel('x(degrees)')
    ylabel('time(ms)')
    zlabel('D_x_t')

    D_x_t_OFF_trans = reshape(D_x_t_OFF_trans, 1,length(y_vect)*length(x_vect)*length(tau_vect));

    %%OFF-CENTER SUSTAINABLE
    disp('OFF sustainable')
    kernel_sign = -1;    %OFF

    %%RGCs PARAMETERS
    %rates from litterature (cf."2_Alpha RGCs types")
    spont_rate_OFF_s = 21;  %Hz
    fig_x = 'OFF_sust_spatial.fig';
    fig_t = 'OFF_sust_temporal.fig';

    %%FIT SPIKING RATES CURVES (fig_x & fig_t)
%     x0 = [0.755370077033390 0.677916947281485 0.0631451612119229 ...
%         0.0819508192686487 0.00735368814964831 0.0104262279522926 ...
%         250.085363778680 94.5940263794188 0.602047179367248 3.72637024055030];
%         %dx=0.1 frame_rate=300

%     x0 = [0.713966982936452,0.515853370652151,0.0558104139198858,...
%         0.120198489362861,0.00715007627820252,0.00731924856237923,...
%         274.276901063888,123.432336143382,0.594126100634334,3.52220193534037];
%         %previous dx=0.7 frame_rate=60
        
%     x0 = [0.712124819578292,0.500004658090514,0.0509603067390575,...
%         0.107462945286736,0.00747131791190198,0.00654957971051557,...
%         276.857528791805,116.791153142591,0.612843798570696,4.06523700365834];
%         %dx=0.7 frame_rate=60

%     x0 = [0.721555369711342,0.500035752189135,0.0566004681420754,...
%         0.155240153522796,0.00852764406858357,0.0147452657317711,...
%         270.042337087629,120.953303137842,0.631324354748957,3.56448863878150];
%         %dx=0.7 frame_rate=120
        
    x0 = [0.733064796158070,0.500052022625911,0.0595691551825236,...
        0.152473627080344,0.00790186493355898,0.0161672569017300,...
        269.759006470483,121.018018712192,0.663782388059425,3.00295286434458];
        %dx=0.7 frame_rate=180
        
    [B_cen, B_sur, alpha_cen, alpha_sur, beta_cen, beta_sur, ...
        A_c, A_s,sigma_c, sigma_s] = ...
        fitSpikingRates(dx,dt,x0,spont_rate_OFF_s,fig_x,fig_t);

    x0 = [B_cen, B_sur, alpha_cen, alpha_sur, beta_cen, beta_sur, ...
        A_c, A_s,sigma_c, sigma_s];

    save('Modelization\x0_OFF_s', 'x0')

    %%DESIGN THE SPATIAL FILTER (kernel D_x)
    D_x_cen_OFF_s = A_c*exp(-(xx.^2 + yy.^2)/...
        (2*(sigma_c^2)))/(2*pi*sigma_c^2);
    D_x_sur_OFF_s = A_s*exp(-(xx.^2 + yy.^2)/...
        (2*(sigma_s^2)))/(2*pi*sigma_s^2);
    D_x = kernel_sign * (D_x_cen_OFF_s - D_x_sur_OFF_s);
    figure
    surf(xx,yy,D_x)
    title('Kernel OFF-sust')
    xlabel('x(degrees)')
    ylabel('y (degrees)')
    zlabel('D_x')

    %%DESIGN TEMPORAL FILTER
    D_t_cen_OFF_s = alpha_cen^2*tau_vect.*exp(-alpha_cen*tau_vect) - ...
       B_cen*beta_cen^2*tau_vect.*exp(-beta_cen*tau_vect);
    D_t_sur_OFF_s = alpha_sur^2*tau_vect.*exp(-alpha_sur*tau_vect) - ...
       B_sur*beta_sur^2*tau_vect.*exp(-beta_sur*tau_vect);
    figure
    subplot(2,1,1)
    plot(tau_vect,D_t_cen_OFF_s)
    title('D-t-cen-OFF-sust')
    xlabel('time(ms)')
    ylabel('D_t')
    subplot(2,1,2)
    plot(tau_vect,D_t_sur_OFF_s)
    title('D-t-sur-OFF-sust')
    xlabel('time(ms)')
    ylabel('D_t')


    %%DESIGN THE SPATIO-TEMPORAL FILTER (kernel D_x_t)
    D_x_t_OFF_sust = zeros(length(y_vect),length(x_vect),length(tau_vect));
    for k = 1:length(tau_vect)
        D_x_t_OFF_sust(:,:,k) = kernel_sign * (D_t_cen_OFF_s(k)*D_x_cen_OFF_s -...
            D_t_sur_OFF_s(k)*D_x_sur_OFF_s);
    end
    [txt,xxt] = meshgrid(tau_vect,x_vect);
    plot_D_x_t = permute(D_x_t_OFF_sust, [2 3 1]);
    figure
    [~,idx_y] = (min(abs(yy(:,1))));
    surf(xxt,txt,plot_D_x_t(:,:,idx_y))
    title('D-x-t-OFF-sust')
    xlabel('x(degrees)')
    ylabel('time(ms)')
    zlabel('D_x_t')

    D_x_t_OFF_sust = reshape(D_x_t_OFF_sust, 1,length(y_vect)*length(x_vect)*length(tau_vect));

    %%


    %start computations
    global temp_step
    temp_step = dt;
%     tau_vect = tau_vect(tau_vect < 500);
%     dt = tau_vect(2)-tau_vect(1);   %converts to ms

    %Downsampling
    pointsPerTrial = (RFmapPx(end)-RFmapPx(1))/dt;
    Step_RFmap = floor(size(RFmapMUA,1)/pointsPerTrial);
    RFmapPx = RFmapPx(1:Step_RFmap:end);
    RFmapMUA = RFmapMUA(1:Step_RFmap:end,:,:);

    RFmapMUA(RFmapPx < -100 | RFmapPx > 500,:,:) = [];
    RFmapPx(RFmapPx < -100 | RFmapPx > 500) = [];

%     RFmapMUA(1,:,:) = [];   %Because the first value of each trial is weird due to filtering
%     RFmapPx(1) = [];

    pointsPerTrial = (RetstimPx(end)-RetstimPx(1))/dt;
    Step_Retstim = floor(size(RetstimMUA,1)/pointsPerTrial);
    RetstimPx = RetstimPx(1:Step_Retstim:end);
    RetstimMUA = RetstimMUA(1:Step_Retstim:end,:,:);
    
    RetstimMUA(RetstimPx < -100,:,:) = [];
    RetstimPx(RetstimPx < -100) = [];

%     RetstimMUA(1,:,:) = [];   %Because the first value of each trial is weird due to filtering
%     RetstimPx(1) = [];

    %% Compute RGCs rates for RFmap
    checksz = 5;    %degree
    small_Mat = zeros(size(RFmapMat,1),3*nb_check);
    small_Mat(:,1:3:end) = RFmapMat(:,2:5:end) + mouseposx;
    small_Mat(:,2:3:end) = RFmapMat(:,3:5:end) + mouseposy;
    small_Mat(:,3:3:3*nb_check/2) = 1;
    small_Mat(:,3*nb_check/2+3:3:end) = -1;

    t_vect = RFmapPx(RFmapPx>=0);

    save_rates = zeros(length(RFmapPx),size(RFmapMUA,2),(2*nb_RGCs+1)*(2*nb_RGCs+1)*4);
    final_integral = zeros(length(t_vect),size(RFmapMUA,2),(2*nb_RGCs+1)*(2*nb_RGCs+1)*4);
    % stim = zeros(size(yy,1),size(xx,2),length(trial),(2*nb_RGCs+1)^2);


    disp('Compute RFmap RGCs rates')
    tic

    for k = 1:size(save_rates,3)/4
        %k=1 => top left | k=2 => 1st column, 2nd row | k=end => bottom right
%         disp(k)
        alphaposx = (ceil(k/(2*nb_RGCs+1))-(nb_RGCs+1));  %from -nb_RGCs to +nb_RGCs
        alphaposy = mod(k,2*nb_RGCs+1);
        if alphaposy == 0
            alphaposy = 2*nb_RGCs+1;  %from 1 to 2*nb_RGCs+1
        end
        alphaposy = - (alphaposy-(nb_RGCs+1));  %from +nb_RGCs to -nb_RGCs
        for trial = 1:size(RFmapMat,1)
            Stim1Trial = zeros(size(yy,1),size(xx,2));
            for i = 0:nb_check-1
                Stim1Trial(...
                    xx>small_Mat(trial,1+3*i)-checksz/2 - alphaposx*alpha_size &...
                    xx<small_Mat(trial,1+3*i)+checksz/2 - alphaposx*alpha_size & ...
                    yy>small_Mat(trial,2+3*i)-checksz/2 - alphaposy*alpha_size &...
                    yy<small_Mat(trial,2+3*i)+checksz/2 - alphaposy*alpha_size )...
                    = small_Mat(trial,3+3*i);
            end
            Stim1Trial(...
                xx > mouseposx+max_x-alphaposx*alpha_size | ...
                xx < mouseposx+min_x-alphaposx*alpha_size | ...
                yy > mouseposy+max_y-alphaposy*alpha_size | ...
                yy < mouseposy+min_y-alphaposy*alpha_size )...
                = 0;
    %         figure
    %         surface(xx,yy,Stim1Trial)
    %         pause()
            for type = 1:4
                switch type
                    case 1
                        kernel_sign = 1;
                        D_x_cen = D_x_cen_ON_t;
                        D_x_sur = D_x_sur_ON_t;
                        D_t_cen = D_t_cen_ON_t;
                        D_t_sur = D_t_sur_ON_t;
                    case 2
                        kernel_sign = 1;
                        D_x_cen = D_x_cen_ON_s;
                        D_x_sur = D_x_sur_ON_s;
                        D_t_cen = D_t_cen_ON_s;
                        D_t_sur = D_t_sur_ON_s;
                    case 3
                        kernel_sign = -1;
                        D_x_cen = D_x_cen_OFF_t;
                        D_x_sur = D_x_sur_OFF_t;
                        D_t_cen = D_t_cen_OFF_t;
                        D_t_sur = D_t_sur_OFF_t;
                    case 4
                        kernel_sign = -1;
                        D_x_cen = D_x_cen_OFF_s;
                        D_x_sur = D_x_sur_OFF_s;
                        D_t_cen = D_t_cen_OFF_s;
                        D_t_sur = D_t_sur_OFF_s;
                end
                spatial_center_integral = sum(sum(dx^2*D_x_cen.*Stim1Trial));
                spatial_surround_integral = sum(sum(dx^2*D_x_sur.*Stim1Trial));

                total_integral = dt*(kernel_sign * ...
                        (D_t_cen*spatial_center_integral - ...
                        D_t_sur*spatial_surround_integral));

                final_integral(1:length(total_integral),trial,...
                    k+(type-1)*size(save_rates,3)/4) = cumsum(total_integral);
                final_integral(length(total_integral)+1:end,trial,...
                    k+(type-1)*size(save_rates,3)/4) = sum(total_integral);
            end
    %             if and(k==25,trial==290)
    %         if and(k==25,trial==62)
    %             spat_cen = spatial_center_integral
    %             spat_sur = spatial_surround_integral
    %             integ = total_integral;
    %             figure
    %             surface(xx,yy,Stim1Trial); colormap gray
    %         end
    %         stim(:,:,trial,k) = Stim1Trial;
        end
    end
    toc

    %% Compute RGCs rates for Retstim
    % Retstim_integral = zeros(length(RetstimPx)-length(RetstimPx(RetstimPx<0)),size(RetstimMUA,2),(2*nb_RGCs+1)*(2*nb_RGCs+1)*4);
    % stim = zeros(size(yy,1),size(xx,2),length(trial),(2*nb_RGCs+1)^2);

    Retstim_rates = zeros(length(RetstimPx),size(RetstimMUA,2),(2*nb_RGCs+1)*(2*nb_RGCs+1)*4);

    t_vect = RetstimPx(RetstimPx>=0);

    disp('Compute Retstim RGCs rates')
    tic

    clear final_stim
    for i = 1:size(Retstim_rates,3)/4
        %k=1 => top left | k=2 => 1st column, 2nd row | k=end => bottom right
%         disp(i)
        alphaposx = (ceil(i/(2*nb_RGCs+1))-(nb_RGCs+1))  %from -nb_RGCs to +nb_RGCs
        alphaposy = mod(i,2*nb_RGCs+1);
        if alphaposy == 0
            alphaposy = 2*nb_RGCs+1;  %from 1 to 2*nb_RGCs+1
        end
        alphaposy = - (alphaposy-(nb_RGCs+1))  %from +nb_RGCs to -nb_RGCs
        period = [1 1];
        my_time = tic;
        for trial = 1:size(RetstimMUA,2)
%             disp('trial')
%             disp(trial)
%             clear RetstimMat
            period(2) = find(strcmp(My_stim.action(period(1):end),'Disappear'), 1) + period(1);
        %     t_vect = My_stim.time(period(1)):dt:(My_stim.time(period(2)) - 1E-4);
            RetstimMat = My_stim(period(1):period(2)-1, :);
            RetstimMat.time = RetstimMat.time - RetstimMat.time(1);
                        

            if i==1
                IDX_STIM = IDX_STIM + 1;
                ALL_STIM(IDX_STIM) = struct('session',session, 'trial',trial, 'stim', My_stim(period(1)+1,:));
            end

            Target_t_vect_long = zeros(size(yy,1),size(xx,2),length(tau_vect));
    %         final_stim = zeros(size(yy,1)*size(xx,2)*length(tau_vect),length(t_vect));

            prev_k = 0;
%                 my_time = tic;
            for t = 1:length(t_vect)
    %             stim = zeros(size(yy,1),size(xx,2));
                k = find(RetstimMat.time <= t_vect(t), 1, 'last');
                if ~(prev_k == k)
                    stim = zeros(size(xx));
                    stim(((xx-RetstimMat.centerx(k) + alphaposx*alpha_size).^2 + ...
                        (yy-RetstimMat.centery(k) + alphaposy*alpha_size).^2) <= ...
                        (RetstimMat.size(k)/2)^2) = RetstimMat.sprite(k);
                    %"(x+xc)^2" instead of "(x-xc)^2" because we
                    %simulate alpha position from stimulus motion
                    stim(...
                        xx > mouseposx+max_x-alphaposx*alpha_size | ...
                        xx < mouseposx+min_x-alphaposx*alpha_size | ...
                        yy > mouseposy+max_y-alphaposy*alpha_size | ...
                        yy < mouseposy+min_y-alphaposy*alpha_size )...
                        = 0; %Limit the stimulus respecting the size of the screen
                    prev_k = k;
    %                 figure
    %                 surface(xx,yy,stim) 
    %                 t
    %                 pause()
    %             else
    %                 disp('same k')
    %                 t
                end
                Target_t_vect_long(:,:,1) = stim;
                final_stim(:,t) = reshape(Target_t_vect_long, size(yy,1)*size(xx,2)*length(tau_vect),1);
                Target_t_vect_long = circshift(Target_t_vect_long,[0 0 1]);
            end
%                 toc(my_time)
            for type = 1:4
                switch type
                    case 1
                        spont_rate = spont_rate_ON_t;
                        D_x_t = D_x_t_ON_trans;
                    case 2
                        spont_rate = spont_rate_ON_s;
                        D_x_t = D_x_t_ON_sust;
                    case 3
                        spont_rate = spont_rate_OFF_t;
                        D_x_t = D_x_t_OFF_trans;
                    case 4
                        spont_rate = spont_rate_OFF_s;
                        D_x_t = D_x_t_OFF_sust;
                end
                spike_rate = repmat(spont_rate,1,length(t_vect)) + dt*dx^2 * (D_x_t*final_stim);
    %                 spike_rate = spont_rate + sum(dt*sum(sum(dx^2*D_x_t.*...
    %                     Target_t_vect_long)))
                %implement the threshold
                spike_rate = max(0,spike_rate);
                %save insight final variable
                Retstim_rates(length(RetstimPx(RetstimPx<0))+1:end,trial,...
                    i+(type-1)*size(Retstim_rates,3)/4) = spike_rate;
    %             pause()
            end  
            
            period(1) = period(2);
        end
        toc(my_time)
    end
    toc
    % Retstim_rates(length(t_vect(t_vect<0))+1:end,:,:) = spont_rate + final_integral;

    %%
    for type = 1:4
        switch type
            case 1
                spont_rate = spont_rate_ON_t;
            case 2
                spont_rate = spont_rate_ON_s;
            case 3
                spont_rate = spont_rate_OFF_t;
            case 4
                spont_rate = spont_rate_OFF_s;
        end
        save_rates(1:length(RFmapPx(RFmapPx<0)),:,...
            1 + (type-1)*size(save_rates,3)/4:type*size(save_rates,3)/4) ...
            = spont_rate;
        save_rates(length(RFmapPx(RFmapPx<0))+1:end,:,...
            1 + (type-1)*size(save_rates,3)/4:type*size(save_rates,3)/4) ...
            = spont_rate + final_integral(:,:,...
            1 + (type-1)*size(save_rates,3)/4:type*size(save_rates,3)/4);

%         Retstim_rates(1:length(RetstimPx(RetstimPx<0)),:,...
%             1 + (type-1)*size(Retstim_rates,3)/4:type*size(Retstim_rates,3)/4)...
%             = spont_rate;
    end
    save_rates = max(0,save_rates);
%     Retstim_rates = max(0,Retstim_rates);
    
    IDX_RGC_RFmap = IDX_RGC_RFmap + 1;
    ALL_RGC_RFmap(IDX_RGC_RFmap) = struct('session',session, 'Rates',save_rates);
%     IDX_RGC_RETSTIM = IDX_RGC_RETSTIM + 1;
%     ALL_RGC_RETSTIM(IDX_RGC_RETSTIM) = struct('session',session, 'Rates',Retstim_rates);   
    
end

save('Results\ALL_STIM.mat','ALL_STIM')
save('Results\ALL_RGC_RFmap.mat','ALL_RGC_RFmap')
save('Results\ALL_RGC_RETSTIM.mat','ALL_RGC_RETSTIM')

%%
%% Compute the SC rates
%% struct to save all results
ALL_WEIGHTS = struct('session',1, 'channel',1, 'gain',1, 't_shift',40, 'spont',40, 'weights',zeros(7));
IDX_WEIGHTS = 0;
ALL_RFmap_RATES = struct('session',1, 'channel',1, 'rates',zeros(500));
IDX_RFmap = 0;
ALL_wrong_WEIGHTS = struct('session',1, 'channel',1, 'gain',1, 't_shift',40, 'spont',40, 'weights',zeros(7));
IDX_wrong_WEIGHTS = 0;
ALL_wrong_RFmap_RATES = struct('session',1, 'channel',1, 'rates',zeros(500));
IDX_wrong_RFmap = 0;
ALL_RETSTIM_RATES = struct('session',1, 'channel',1, 'rates',zeros(500));
IDX_RETSTIM = 0;

%%
for session = 1:5   %5 different locations recorded
    disp('session')
    disp(session)
    close all
    
%     dx = .7;    %degree
%     frame_rate = 60*3;   %Hz
%     dt = 1/frame_rate * 1000; %ms
%     global temp_step
%     temp_step = dt;

    %% Load data
    %%Load RFmap computed rates
    my_path = 'Results\ALL_RGC_RFmap.mat';
    ALL_RGC_RFmap = cell2mat(struct2cell(load(my_path)));
    save_rates = ALL_RGC_RFmap(session).Rates;
    %%Load Retstim stimuli
    my_path = 'Results\ALL_STIM.mat';
    ALL_STIM = cell2mat(struct2cell(load(my_path)));
%     Retstim_rates = ALL_RGC_RETSTIM(session).Rates;
    %%Load Retstim computed rates
    my_path = 'Results\ALL_RGC_RETSTIM.mat';
    ALL_RGC_RETSTIM = cell2mat(struct2cell(load(my_path)));
    Retstim_rates = ALL_RGC_RETSTIM(session).Rates;
    %%Load RFmap log
    my_path = ['Recordings\Session-', int2str(session), '\RFmapMat.mat'];
    RFmapMat = cell2mat(struct2cell(load(my_path)));
    %%Load Multi-chanel spike rates (RFmap)
    my_path = ['Recordings\Session-', int2str(session), '\RFmapMUA.mat'];
    RFmapMUA = cell2mat(struct2cell(load(my_path)));
    %%Load time vector (RFmap)
    my_path = ['Recordings\Session-', int2str(session), '\RFmapPx.mat'];
    RFmapPx = cell2mat(struct2cell(load(my_path)));
    %%Load details (RFmap)
    my_path = ['Recordings\Session-', int2str(session), '\details.mat'];
    details = cell2mat(struct2cell(load(my_path)));
    %%Load Multi-chanel spike rates (Retstim)
    my_path = ['Recordings\Session-', int2str(session), '\RetstimMUA.mat'];
    RetstimMUA = cell2mat(struct2cell(load(my_path)));
    %%Load Retstim px (time vector)
    my_path = ['Recordings\Session-', int2str(session), '\RetStimPx.mat'];
    RetstimPx = cell2mat(struct2cell(load(my_path)));
    
    %% Average the RFmapMUA rates (actual and computed)
    switch session
        case 1
            %session 1 x=55,y=15
            x_ref = 55-10:5:55+10;
            y_ref = 15-10:5:15+10;
        case 2
            %session 2 x=15,y=25
            x_ref = 15-10:5:15+10;
            y_ref = 25-10:5:25+10;
        case 3
            %session 3 x=15,y=10
            x_ref = 15-10:5:15+10;
            y_ref = 10-10:5:10+10;
        case 4
            %session 4 x=-30, y=45
            x_ref = -30-10:5:-30+10;
            y_ref = 45-10:5:45+10;
        case 5
            %session 5 x=30, y=40
            x_ref = 30-10:5:30+10;
            y_ref = 40-10:5:40+10;
    end
    
    clear mean_RFmapMUA
    clear mean_RGC_RFmap

    idx_mean = 0;
    for color = 1:2
        for x = x_ref
            for y = y_ref
                [row,~] = find(details(:,2)==x & details(:,3)==y);
                if ~isempty(row)
                    ID = details(row,1);
                    [row, col] = find(and(RFmapMat(:,1:5:60) == ID, RFmapMat(:,5:5:60) == color));
                    row=sort(row);
                    if ~isempty(row)
                        idx_mean = idx_mean + 1;
                        mean_RFmapMUA(:,idx_mean,:) = mean(RFmapMUA(:,row,:),2);
                        mean_RGC_RFmap(:,idx_mean,:) = mean(save_rates(:,row,:),2);
                        if idx_mean==10
                            disp(ID)
                        end

                    end
                end
            end
        end
    end
    
    %% Smooth
    my_smooth = smooth(mean_RFmapMUA(:),30);
    smooth_RFmapMUA = reshape(my_smooth,size(mean_RFmapMUA,1),...
        size(mean_RFmapMUA,2),size(mean_RFmapMUA,3));

    %% Downsampling
    pointsPerTrial = (RFmapPx(end)-RFmapPx(1))/dt;
    Step_RFmap = floor(size(smooth_RFmapMUA,1)/pointsPerTrial);
    RFmapMUA = smooth_RFmapMUA(1:Step_RFmap:end,:,:);
    RFmapPx = RFmapPx(1:Step_RFmap:end);
    
    RFmapMUA(RFmapPx < -100 | RFmapPx > 500,:,:) = [];
    RFmapPx(RFmapPx < -100 | RFmapPx > 500) = [];

    

    pointsPerTrial = (RetstimPx(end)-RetstimPx(1))/dt;
    Step_Retstim = floor(size(RetstimMUA,1)/pointsPerTrial);
    RetstimMUA = RetstimMUA(1:Step_Retstim:end,:,:);
    RetstimPx = RetstimPx(1:Step_Retstim:end);
    
    RetstimMUA(RetstimPx < -100,:,:) = [];
    RetstimPx(RetstimPx < -100) = [];
    
    %% Wrong mean_RGC_RFmap (by reversing the time axis) for the comparison
    wrong_RFmap_rates = flipud(mean_RGC_RFmap);
    
    %% Start computations
    My_chan = [4:7 17:21];  %channels in superficial(17:21) and in deep(4:7)
%     My_chan = [19 17];  %channels in superficial(17:21) and in deep(4:7)
    for channel = My_chan
        %% Find the optimal weights of each alpha RFG from the RFmap stimuli (~white noise)
        %Assumption of uniformly distributed weights along the "radius"
        %Example: x0 = [.1 .2 .3 .4 .5 .6 .7 .8 .9 1];
        %So, weights = 
        %.1 .2 .4 .7 .4 .2 .1 
        %.2 .3 .5 .8 .5 .3 .2 
        %.4 .5 .6 .9 .6 .5 .4
        %.7 .8 .9  1 .9 .8 .7
        %.4 .5 .6 .9 .6 .5 .4
        %.2 .3 .5 .8 .5 .3 .2
        %.1 .2 .4 .7 .4 .2 .1
        total_time = tic;
        [gain, t_shift, thresh, weights] = computeWeights(mean_RGC_RFmap,RFmapMUA(:,:,channel));
        shift_step = round(t_shift/dt);
        IDX_WEIGHTS = IDX_WEIGHTS + 1;
        ALL_WEIGHTS(IDX_WEIGHTS) = struct('session',session, ...
            'channel',channel, 'gain',gain, 't_shift',t_shift, 'thresh',thresh, 'weights',weights);
        
        [wrong_gain, wrong_t_shift, wrong_thresh, wrong_weights] = computeWeights(wrong_RFmap_rates,RFmapMUA(:,:,channel));
        wrong_shift_step = round(wrong_t_shift/dt);
        IDX_wrong_WEIGHTS = IDX_WEIGHTS + 1;
        ALL_wrong_WEIGHTS(IDX_WEIGHTS) = struct('session',session, ...
            'channel',channel, 'gain',wrong_gain, 't_shift',wrong_t_shift, 'thresh',wrong_thresh, 'weights',wrong_weights);
        toc(total_time)

        %% Final rates for RFmap
        shifted_rates = circshift(mean_RGC_RFmap,[shift_step 0 0]);
        shifted_rates(1:shift_step+1,:,:) = repmat(mean_RGC_RFmap(1,:,:),shift_step+1,1,1);
        new_rates = reshape(shifted_rates,size(mean_RGC_RFmap,1)*size(mean_RGC_RFmap,2), size(mean_RGC_RFmap,3)/4, 4);
        new_weights = reshape(weights,size(weights,1)*size(weights,2),size(weights,3));
        weighted_rates = zeros(size(new_rates,1),1);
        for i = 1:4
            weighted_rates = weighted_rates + new_rates(:,:,i) * new_weights(:,i);
        end
        RFmap_final_rates = reshape(-thresh + gain*weighted_rates,size(mean_RGC_RFmap,1),size(mean_RGC_RFmap,2));
        RFmap_final_rates = max(RFmap_final_rates,0);
        IDX_RFmap = IDX_RFmap + 1;
        ALL_RFmap_RATES(IDX_RFmap) = struct('session',session, 'channel',channel, 'rates',RFmap_final_rates);

        %% Final rates for wrong RFmap
        shifted_rates = circshift(wrong_RFmap_rates,[wrong_shift_step 0 0]);
        shifted_rates(1:wrong_shift_step+1,:,:) = repmat(wrong_RFmap_rates(1,:,:),wrong_shift_step+1,1,1);
        new_rates = reshape(shifted_rates,size(wrong_RFmap_rates,1)*size(wrong_RFmap_rates,2), size(wrong_RFmap_rates,3)/4, 4);
        new_weights = reshape(wrong_weights,size(wrong_weights,1)*size(wrong_weights,2),size(wrong_weights,3));
        weighted_rates = zeros(size(new_rates,1),1);
        for i = 1:4
            weighted_rates = weighted_rates + new_rates(:,:,i) * new_weights(:,i);
        end
        RFmap_final_rates = reshape(-wrong_thresh + wrong_gain*weighted_rates,size(wrong_RFmap_rates,1),size(wrong_RFmap_rates,2));
        RFmap_final_rates = max(RFmap_final_rates,0);
        IDX_wrong_RFmap = IDX_RFmap + 1;
        ALL_wrong_RFmap_RATES(IDX_RFmap) = struct('session',session, 'channel',channel, 'rates',RFmap_final_rates);
        
        %% Final rates for Retstim  %%doesn't work for now...
        shifted_rates = circshift(Retstim_rates,shift_step);
        shifted_rates(1:shift_step+1,:,:) = repmat(Retstim_rates(1,:,:),shift_step+1,1,1);
        new_rates = reshape(shifted_rates,size(Retstim_rates,1)*size(Retstim_rates,2), size(Retstim_rates,3)/4, 4);
        new_weights = reshape(weights,size(weights,1)*size(weights,2),size(weights,3));
        weighted_rates = zeros(size(new_rates,1),1);
        for i = 1:4
            weighted_rates = weighted_rates + new_rates(:,:,i) * new_weights(:,i);
        end
        Retstim_final_rates = reshape(-thresh + gain*weighted_rates,size(Retstim_rates,1),size(Retstim_rates,2));
        Retstim_final_rates = max(Retstim_final_rates,0);
        IDX_RETSTIM = IDX_RETSTIM + 1;
        ALL_RETSTIM_RATES(IDX_RETSTIM) = struct('session',session, 'channel',channel, 'rates',Retstim_final_rates);
        
        disp('done')
        disp(channel)
        
    end
    
    % try to put it outside the loop to save everything
    save('Results\REF_RFmapMUA.mat','RFmapMUA')
    save('Results\REF_RestimMUA.mat','RetstimMUA')
    
end


%% Save the final variables
save('Results\ALL_WEIGHTS.mat','ALL_WEIGHTS')
save('Results\ALL_RFmap_RATES.mat','ALL_RFmap_RATES')
save('Results\ALL_RETSTIM_RATES.mat','ALL_RETSTIM_RATES')

















%% Previous
% %% Compute the SC rates
% for session = 1:2   %5 different locations recorded
%     disp('session')
%     disp(session)
%     close all
% 
%     %% Load data
%     %%Load RFmap computed rates
%     my_path = 'Results\ALL_RGC_RFmap.mat';
%     ALL_RGC_RFmap = (load(my_path));
%     save_rates = ALL_RGC_RFmap(session).Rates;
%     %%Load Retstim computed rates
%     my_path = 'Results\ALL_RGC_RETSTIM.mat';
%     ALL_RGC_RETSTIM = (load(my_path));
%     Retstim_rates = ALL_RGC_RETSTIM(session).Rates;
%     %%Load Multi-chanel spike rates (RFmap)
%     my_path = ['Recordings\Session-', int2str(session), '\RFmapMUA.mat'];
%     RFmapMUA = cell2mat(struct2cell(load(my_path)));
%     %%Load details (RFmap)
%     my_path = ['Recordings\Session-', int2str(session), '\details.mat'];
%     details = cell2mat(struct2cell(load(my_path)));
%     %%Load Multi-chanel spike rates (Retstim)
%     my_path = ['Recordings\Session-', int2str(session), '\RetstimMUA.mat'];
%     RetstimMUA = cell2mat(struct2cell(load(my_path)));
%     
%     %% Average the RFmapMUA rates (actual and computed)
%     switch session
%         case 1
%             %session 1 x=55,y=15
%             x_ref = 55-10:5:55+10;
%             y_ref = 15-10:5:15+10;
%         case 2
%             %session 2 x=15,y=25
%             x_ref = 15-10:5:15+10;
%             y_ref = 25-10:5:25+10;
%         case 3
%             %session 3 x=15,y=10
%             x_ref = 15-10:5:15+10;
%             y_ref = 10-10:5:10+10;
%         case 4
%             %session 4 x=-30, y=45
%             x_ref = -30-10:5:-30+10;
%             y_ref = 45-10:5:45+10;
%         case 5
%             %session 5 x=30, y=40
%             x_ref = 30-10:5:30+10;
%             y_ref = 40-10:5:40+10;
%     end
% 
%     idx_mean = 0;
%     for color = 1:2
%         for x = x_ref
%             for y = y_ref
%                 [row,~] = find(details(:,2)==x & details(:,3)==y);
%                 if ~isempty(row)
%                     ID = details(row,1);
%                     [row, col] = find(and(RFmapMat(:,1:5:60) == ID, RFmapMat(:,5:5:60) == color));
%                     row=sort(row);
%                     if ~isempty(row)
%                         idx_mean = idx_mean + 1;
%                         mean_RFmapMUA(:,idx_mean,:) = mean(RFmapMUA(:,row,:),2);
%                         mean_RGC_RFmap(:,idx_mean,:) = mean(save_rates(:,row,:),2);
% 
%                     end
%                 end
%             end
%         end
%     end
%     
%     %% Smooth
%     my_smooth = smooth(mean_RFmapMUA(:),30);
%     smooth_RFmapMUA = reshape(my_smooth,size(mean_RFmapMUA,1),...
%         size(mean_RFmapMUA,2),size(mean_RFmapMUA,3));
%     
%     %% Downsampling
%     pointsPerTrial = (RFmapPx(end)-RFmapPx(1))/dt;
%     Step_RFmap = floor(size(smooth_RFmapMUA,1)/pointsPerTrial);
%     RFmapMUA = smooth_RFmapMUA(1:Step_RFmap:end,:,:);
% 
%     RFmapMUA(RFmapPx < -100 | RFmapPx > 500,:,:) = [];
%     
% 
%     pointsPerTrial = (RetstimPx(end)-RetstimPx(1))/dt;
%     Step_Retstim = floor(size(RetstimMUA,1)/pointsPerTrial);
%     RetstimMUA = RetstimMUA(1:Step_Retstim:end,:,:);
%     
%     RetstimMUA(RetstimPx < -100) = [];
%     
%     %% Wrong mean_RGC_RFmap (by reversing the time axis) for the comparison
%     wrong_RFmap_rates = flipud(mean_RGC_RFmap);
%     
%     
% 
%     %% Start computations
% %     My_chan = [4:7 17:21];  %channels in superficial(17:21) and in deep(4:7)
%     My_chan = [7 17];  %channels in superficial(17:21) and in deep(4:7)
%     for channel = My_chan
%         %% Find the optimal weights of each alpha RFG from the RFmap stimuli (~white noise)
%         %Assumption of uniformly distributed weights along the "radius"
%         %Example: x0 = [.1 .2 .3 .4 .5 .6 .7 .8 .9 1];
%         %So, weights = 
%         %.1 .2 .4 .7 .4 .2 .1 
%         %.2 .3 .5 .8 .5 .3 .2 
%         %.4 .5 .6 .9 .6 .5 .4
%         %.7 .8 .9  1 .9 .8 .7
%         %.4 .5 .6 .9 .6 .5 .4
%         %.2 .3 .5 .8 .5 .3 .2
%         %.1 .2 .4 .7 .4 .2 .1
%         total_time = tic;
%         [gain, t_shift, spont, weights] = computeWeights(mean_RGC_RFmap,RFmapMUA(:,:,channel));
%         shift_step = round(t_shift/dt);
%         IDX_WEIGHTS = IDX_WEIGHTS + 1;
%         ALL_WEIGHTS(IDX_WEIGHTS) = struct('session',session, ...
%             'channel',channel, 'gain',gain, 't_shift',t_shift, 'spont',spont, 'weights',weights);
%         
%         [wrong_gain, wrong_t_shift, wrong_spont, wrong_weights] = computeWeights(wrong_RFmap_rates,RFmapMUA(:,:,channel));
%         wrong_shift_step = round(wrong_t_shift/dt);
%         IDX_wrong_WEIGHTS = IDX_WEIGHTS + 1;
%         ALL_wrong_WEIGHTS(IDX_WEIGHTS) = struct('session',session, ...
%             'channel',channel, 'gain',wrong_gain, 't_shift',wrong_t_shift, 'spont',wrong_spont, 'weights',wrong_weights);
%         toc(total_time)
% 
%         %% Final rates for RFmap
%         shifted_rates = circshift(mean_RGC_RFmap,[shift_step 0 0]);
%         shifted_rates(1:shift_step+1,:,:) = repmat(mean_RGC_RFmap(1,:,:),shift_step+1,1,1);
%         new_rates = reshape(shifted_rates,size(mean_RGC_RFmap,1)*size(mean_RGC_RFmap,2), size(mean_RGC_RFmap,3)/4, 4);
%         new_weights = reshape(weights,size(weights,1)*size(weights,2),size(weights,3));
%         weighted_rates = zeros(size(new_rates,1),1);
%         for i = 1:4
%             weighted_rates = weighted_rates + new_rates(:,:,i) * new_weights(:,i);
%         end
%         RFmap_final_rates = reshape(spont + gain*weighted_rates,size(mean_RGC_RFmap,1),size(mean_RGC_RFmap,2));
%         IDX_RFmap = IDX_RFmap + 1;
%         ALL_RFmap_RATES(IDX_RFmap) = struct('session',session, 'channel',session, 'rates',RFmap_final_rates);
% 
%         %% Final rates for wrong RFmap
%         shifted_rates = circshift(wrong_RFmap_rates,[wrong_shift_step 0 0]);
%         shifted_rates(1:wrong_shift_step+1,:,:) = repmat(wrong_RFmap_rates(1,:,:),wrong_shift_step+1,1,1);
%         new_rates = reshape(shifted_rates,size(wrong_RFmap_rates,1)*size(wrong_RFmap_rates,2), size(wrong_RFmap_rates,3)/4, 4);
%         new_weights = reshape(wrong_weights,size(wrong_weights,1)*size(wrong_weights,2),size(wrong_weights,3));
%         weighted_rates = zeros(size(new_rates,1),1);
%         for i = 1:4
%             weighted_rates = weighted_rates + new_rates(:,:,i) * new_weights(:,i);
%         end
%         RFmap_final_rates = reshape(wrong_spont + wrong_gain*weighted_rates,size(wrong_RFmap_rates,1),size(wrong_RFmap_rates,2));
%         IDX_wrong_RFmap = IDX_RFmap + 1;
%         ALL_wrong_RFmap_RATES(IDX_RFmap) = struct('session',session, 'channel',session, 'rates',RFmap_final_rates);
%         
%         %% Final rates for Retstim
%         shifted_rates = circshift(Retstim_rates,shift_step);
%         shifted_rates(1:shift_step+1,:,:) = repmat(Retstim_rates(1,:,:),shift_step+1,1,1);
%         new_rates = reshape(shifted_rates,size(Retstim_rates,1)*size(Retstim_rates,2), size(Retstim_rates,3)/4, 4);
%         new_weights = reshape(weights,size(weights,1)*size(weights,2),size(weights,3));
%         weighted_rates = zeros(size(new_rates,1),1);
%         for i = 1:4
%             weighted_rates = weighted_rates + new_rates(:,:,i) * new_weights(:,i);
%         end
%         Retstim_final_rates = reshape(spont + gain*weighted_rates,size(Retstim_rates,1),size(Retstim_rates,2));
%         IDX_RETSTIM = IDX_RETSTIM + 1;
%         ALL_RETSTIM_RATES(IDX_RETSTIM) = struct('session',session, 'channel',channel, 'rates',Retstim_final_rates);
%         
%         disp('done')
%         disp(channel)
%         
%     end
%     
%     
% end
% 
% 
% %% Save the final variables
% save('Results\REF_RFmapMUA.mat','RFmapMUA')
% save('Results\REF_RestimMUA.mat','RetstimMUA')
% save('Results\ALL_WEIGHTS.mat','ALL_WEIGHTS')
% save('Results\ALL_RFmap_RATES.mat','ALL_RFmap_RATES')
% save('Results\ALL_RETSTIM_RATES.mat','ALL_RETSTIM_RATES')
% 
% 
% 



%%
%     final_rates = max(0,final_rates);

%     figure
%     plot(t_vect,final_rates)
%     title('Computed final rates')
%     xlabel('time(ms)')
%     ylabel('Spike rate(Hz)')


    % x0 = (1:sum(1:ceil(length(weights)/2)))/sum(1:ceil(length(weights)/2));
%     x0 = rand(sum(1:ceil(length(weights)/2)),1);

%     new_weights = computeWeights(t_vect,save_rates,final_rates);


%%
%%



% t_vect = RetstimPx;
% save_rates = zeros(length(t_vect),size(RetstimMUA,2),(2*nb_RGCs+1)*(2*nb_RGCs+1)*4);
% save_rates(1:length(t_vect(t_vect<0)),:,:) = spont_rate;
% RetstimMat = zeros(length(t_vect),size(My_stim,2),size(RetstimMUA,2));
% 
% trial = 0;
% period = [1 1];
% while (period(2) < height(My_stim))
%     period(2) = find(strcmp(My_stim.action(period(1):end),'Disappear'), 1) + period(1);
% %     t_vect = My_stim.time(period(1)):dt:(My_stim.time(period(2)) - 1E-4);
%     period(1) = period(2);
%     
% 
% %     save_rates = zeros(2*nb_RGCs+1,2*nb_RGCs+1,length(t_vect));
%     trial = trial+1;
%     i=0;
%     for k1 = -nb_RGCs:nb_RGCs
%         i = i+1;
%         j = 0;
%         for k2 = -nb_RGCs : nb_RGCs
%             j = j+1;
% 
% 
% %% TO BE DONE
% %     for k = 1:size(save_rates,3)/4
% %         %k=1 => top left | k=2 => 1st column, 2nd row | k=end => bottom right
% %         disp(k)
% %         alphaposx = (ceil(k/(2*nb_RGCs+1))-(nb_RGCs+1))  %from -nb_RGCs to +nb_RGCs
% %         alphaposy = mod(k,2*nb_RGCs+1);
% %         if alphaposy == 0
% %             alphaposy = 2*nb_RGCs+1;  %from 1 to 2*nb_RGCs+1
% %         end
% %         alphaposy = - (alphaposy-(nb_RGCs+1))  %from +nb_RGCs to -nb_RGCs
% 
% %%
%     
%     
%     
%             Target_t_vect_long = zeros(size(xx,1),size(xx,1),length(tau_vect));
%             stim = zeros(size(xx));
%             prev_k = 0;
%             for t = 1:length(t_vect)
%                 Target_t_vect_long = circshift(Target_t_vect_long,[0 0 1]);
%                 k = find(My_stim.time <= t_vect(t), 1, 'last');
%                 if ~(My_stim.sprite(k) == 0)
% 
%                     stim = zeros(size(xx));
%                     stim(((xx-My_stim.centerx(k) + k1*alpha_size).^2 + ...
%                         (yy-My_stim.centery(k) + k2*alpha_size).^2) <= ...
%                         (My_stim.size(k)/2)^2) = My_stim.sprite(k);
%                     %"(x+xc)^2" instead of "(x-xc)^2" because we
%                     %simulate alpha position from stimulus motion
%                     stim(xx > screenx/2-stim_centerx - k1*alpha_size | ...
%                         xx < -screenx/2-stim_centerx - k1*alpha_size | ...
%                         yy > screeny/2-stim_centery - k2*alpha_size | ...
%                         yy < -screeny/2-stim_centery - k2*alpha_size) = 0;
%                     %Limit the stimulus respecting the size of the screen
% %                     if ~(k2 == prev_k)
% %                         figure
% %                         surface(xx,yy,stim) 
% %                         t
% %                         pause()
% %                         prev_k = k2;
% %                     end
%                 end
%                 Target_t_vect_long(:,:,1) = stim;
%                 spike_rate = spont_rate + sum(dt*sum(sum(dx^2*D_x_t_ON_trans.*...
%                     Target_t_vect_long)));
%                 %implement the threshold
%                 spike_rate = max(0,spike_rate);
%                 save_rates(i,j,t) = spike_rate;
%             end
%         end
%     end
    
%     weights_size = size(save_rates);
%     weights = zeros(weights_size(1:2));
%     i=0;
%     for k1 = -nb_RGCs:nb_RGCs
%         i = i+1;
%         j = 0;
%         for k2 = -nb_RGCs : nb_RGCs
%             j = j+1;
%             if ~and(k1 == 0, k2 == 0)
%                 weights(i,j) = 1/(k1^2+k2^2);
%                 if(weights(i,j) == 1)
%                     weights(i,j) = 0.9;
%                 end
%             else
%                 weights(i,j)=1;
%             end
%         end
%     end
    % !!! TO BE DONE !!!
% end


%%
% %Compute the rates for each RGC
% t = 0;
% t_vect = px;
% 
% save_rates = zeros(length(t_vect),size(Mua,2),(2*nb_RGCs+1)*(2*nb_RGCs+1));
% save_rates(1:length(t_vect(t_vect<0)),:,:) = spont_rate;
% final_integral = zeros(length(t_vect)-length(t_vect(t_vect<0)),size(Mua,2),(2*nb_RGCs+1)*(2*nb_RGCs+1));
% 
% % i=0;
% % for k1 = -nb_RGCs:nb_RGCs
% %     i = i+1;
% %     j = 0;
% %     for k2 = -nb_RGCs : nb_RGCs
% %         j = j+1;
% disp('begin compute save_rates')
% tic
% for k = 1:size(stim,4)
%     for trial = 1:size(Mat,1)
%         spatial_center_integral = sum(sum(dx^2*D_x_cen.*stim(:,:,trial,k)));
%         spatial_surround_integral = sum(sum(dx^2*D_x_sur.*stim(:,:,trial,k)));
% 
%         total_integral = dt*(kernel_sign * ...
%                 (D_t_cen*spatial_center_integral - ...
%                 D_t_sur*spatial_surround_integral));
% 
%         final_integral(1:length(total_integral),trial,k) = cumsum(total_integral);
% %         if length(total_integral) < (length(t_vect)-length(t_vect(t_vect<0)))
%         final_integral(length(total_integral)+1:end,trial,k) = sum(total_integral);
% %         end
%         if and(k==25,trial==290)
%             spat_cen = spatial_center_integral
%             pat_sur = spatial_surround_integral
%             integ = total_integral;
%         end
%     end
% end
% toc
% save_rates(length(t_vect(t_vect<0))+1:end,:,:) = spont_rate + final_integral;
% save_rates = max(0,save_rates);

%Find the optimized weights for each alpha RGC

%% testAlphaCombination
%Limit the computations at the size of the screen
% x_vect = -(screenx/2+stim_centerx):dx:(screenx/2-stim_centerx);
% y_vect = -(screeny/2+stim_centery):dx:(screeny/2-stim_centery);
% [xx, yy] = meshgrid(x_vect,y_vect);
% yy = flipud(yy);    %To keep the same plane of coordinates

% My_stim.centery(1) = screeny/2-stim_centery;