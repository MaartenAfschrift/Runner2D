classdef OptSpeedRunner_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        SimulateButton                matlab.ui.control.Button
        VideoButton                   matlab.ui.control.Button
        ResultsTableButton            matlab.ui.control.Button
        SimulatonoutputTextAreaLabel  matlab.ui.control.Label
        SimulationResults             matlab.ui.control.TextArea
        FirstNameEditFieldLabel       matlab.ui.control.Label
        FirstNameEditField            matlab.ui.control.EditField
        DagvandewetenschapLabel       matlab.ui.control.Label
        Image                         matlab.ui.control.Image
        ClearWindowButton             matlab.ui.control.Button
        ScoreAxis                     matlab.ui.control.UIAxes
        TendonstiffnessPanel          matlab.ui.container.Panel
        NormachillestendonstiffnessSliderLabel  matlab.ui.control.Label
        NormachillestendonstiffnessSlider  matlab.ui.control.Slider
        UIAxes                        matlab.ui.control.UIAxes
        Image3                        matlab.ui.control.Image
        MuscleforcePanel              matlab.ui.container.Panel
        Button                        matlab.ui.control.Button
        Button_2                      matlab.ui.control.Button
        Button_3                      matlab.ui.control.Button
        Button_4                      matlab.ui.control.Button
        Button_5                      matlab.ui.control.Button
        Button_6                      matlab.ui.control.Button
        MuscleForceNEditFieldLabel    matlab.ui.control.Label
        MuscleForceNEditField         matlab.ui.control.NumericEditField
        SoleusEditFieldLabel          matlab.ui.control.Label
        SoleusEditField               matlab.ui.control.NumericEditField
        GastrocnemiusEditFieldLabel   matlab.ui.control.Label
        GastrocnemiusEditField        matlab.ui.control.NumericEditField
        RectusfemorisEditFieldLabel   matlab.ui.control.Label
        RectusfemorisEditField        matlab.ui.control.NumericEditField
        MuscleForceLabel              matlab.ui.control.Label
        Image2                        matlab.ui.control.Image
        ReadyLampLabel                matlab.ui.control.Label
        ReadyLamp                     matlab.ui.control.Lamp
        MusculoskeletalmodelLabel     matlab.ui.control.Label
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            global MainPath
            global S
            global MForceAv
            global vis
            global osimModel
            global model_state
            
            % add path to simulation code
            TxtPath = which('MyID_PredSim2D_ab87ex.txt');
            [MainPath,~,~] = fileparts(TxtPath);
            addpath(genpath(fullfile(MainPath)));
            app.SimulationResults.Value = {'Start GUI'};
            
            % add the casadipath
            CasFilePath = which('casadiMEX.mexw64');
            [CasPath,~,~] = fileparts(CasFilePath);
            addpath(genpath(CasPath));
            
            % default simulation settings
            MForceAv = app.MuscleForceNEditField.Value;
            
            % something about building model here (from global structure
            % with model input information)
            % Run the simulations
            S.v_tgt     = 1.3;        % average speed
            S.N         = 50;       % number of mesh intervals
            
            % path repository
            S.pathRepo = MainPath;
            
            % quasi random initial guess, pelvis y position
            S.IG_PelvisY = 0.93;   % subject 1 poggensee
            
            % output folder
            S.ResultsFolder     = 'DagWetenschap';
            
            % selection folder with Casadi Functions
            S.CasadiFunc_Folders = 'Gait18_shorterHamstring';
            
            % select folder with polynomials
            S.PolyFolder ='Gait18_shorterHamstring';
            
            % initial guess based on simulations without exoskeletons
            S.IDig          = 2;        % (1) == quasi random, (2) == IK fisdle, (3) from solution
            S.ResultsF_ig   = 'DagWetenschap';
            S.savename_ig   = 'IG_DatWetenschap2';
            
            % Select model
            S.ModelName = 'Gait18'; % other option is 'Rajagopal'
            
            % mass of the subject
            S.mass = 64;    % mass in kg
            
            % default stiffness of mtp joint
            S.k_mtp = 100;
            
            % linear solver
            S.linear_solver = 'mumps';
            
            % Simulation without exoskeelton
            %S.ExternalFunc  = 'FastRunner.dll';
            S.ExternalFunc  = 'PredSim_2D_default_mtp_v3.dll';
            S.savename      = 'DefaultSim';
            
            % tolerance ipopt
            S.tol_ipopt   = 5;    % tolerance ipopt
            
            % mot file with results
            S.MotFile = fullfile(S.pathRepo,'Results',S.ResultsFolder,[S.savename '_q.mot']); % file with outputs
            
            % initialise the figure
            NormachillestendonstiffnessSliderValueChanged(app,true);
            
            % initialise the opensim visualisation
            import org.opensim.modeling.*
            modelPath = fullfile(S.pathRepo,'App','Gait18_Visual.osim');
            osimModel = Model(modelPath);
            osimModel.setUseVisualizer(true);
            model_state = osimModel.initSystem();
            vis = osimModel.getVisualizer();
            vis.addDirToGeometrySearchPaths(char(fullfile(S.pathRepo,'App','Geometry')));
            pause(0.1);
            vis.show(model_state);
        end

        % Button pushed function: SimulateButton
        function SimulateButtonPushed(app, event)
            % callback of simulate function
            %global MainPath
            global Results
            global S
            persistent ScoreVect
            
            % adapt model based on user input
            app.SimulationResults.Value = [app.SimulationResults.Value; {'Adapting strength of the muscles in the model'}];
            
            % import default model
            S.ModelPath = fullfile(S.pathRepo,'OpenSimModel','Gait18_shorterHamstring.osim');
            
            % adapt model based on user input
            import org.opensim.modeling.*
            m = Model(S.ModelPath);
            
            % adapt soleus force
            mSel = m.getMuscles.get('soleus_r');
            FSol = mSel.getMaxIsometricForce();
            mSel.setMaxIsometricForce(FSol + app.SoleusEditField.Value);
            mSel = m.getMuscles.get('soleus_l');
            FSol = mSel.getMaxIsometricForce();
            mSel.setMaxIsometricForce(FSol + app.SoleusEditField.Value);
            
            % adapt gastrocnemius force
            mSel = m.getMuscles.get('gastroc_r');
            FSel = mSel.getMaxIsometricForce();
            mSel.setMaxIsometricForce(FSel + app.GastrocnemiusEditField.Value);
            mSel = m.getMuscles.get('gastroc_l');
            FSel = mSel.getMaxIsometricForce();
            mSel.setMaxIsometricForce(FSel + app.GastrocnemiusEditField.Value);
            
            % adapt rectus femoris force
            mSel = m.getMuscles.get('rect_fem_r');
            FSel = mSel.getMaxIsometricForce();
            mSel.setMaxIsometricForce(FSel + app.RectusfemorisEditField.Value);
            mSel = m.getMuscles.get('rect_fem_l');
            FSel = mSel.getMaxIsometricForce();
            mSel.setMaxIsometricForce(FSel + app.RectusfemorisEditField.Value);
            m.print(fullfile(S.pathRepo,'OpenSimModel','ModelDagWetenschap.osim'));
            S.ModelPath = fullfile(S.pathRepo,'OpenSimModel','ModelDagWetenschap.osim');
            
            % create the casadifunctions for this model
            S.CasadiFunc_Folders = 'DagWetenschapTemp';
            SettingsCasFunc.kTendon_CalfM = app.NormachillestendonstiffnessSlider.Value;
            S.ModelName = 'Gait18';
            CreateCasadiFunctions(S.pathRepo, S.ModelName, S.ModelPath, S.CasadiFunc_Folders,S.PolyFolder,SettingsCasFunc,true);
            
            % print start simulation
            app.SimulationResults.Value = [app.SimulationResults.Value; {'Start simulation'}];
            
            % adjust the savename if this file already exists
            BoolNameAdjusted = false;
            if exist(fullfile(S.pathRepo,'Results',S.ResultsFolder,[S.savename '.mat']),'file')
                BoolNameAdjusted = true;
                NameOr = S.savename;
                NameNew = S.savename;
                ct = 1;
                while exist(fullfile(S.pathRepo,'Results',S.ResultsFolder,[NameNew '.mat']),'file')
                    NameNew = [S.savename '_' num2str(ct)];
                    ct = ct+1;
                end
                S.savename = NameNew;
                app.SimulationResults.Value = [app.SimulationResults.Value; {[NameOr ' was already used, file will be saved as ' S.savename]}];
            end
            S.MotFile = fullfile(S.pathRepo,'Results',S.ResultsFolder,[S.savename '_q.mot']); % file with outputs
            app.ReadyLamp.Color = 'red';
            drawnow
            % predictive simulation
            [Results] = f_PredSim_2D_trapezoidal_OptSpeed(S);
            app.ReadyLamp.Color = 'green';
            drawnow
            % adjust savename to original name (this causes problems with
            % the visulaisation
            if BoolNameAdjusted
                S.savename= NameOr;
            end
            if Results.stats.success  % only store result when simulation was succesfull
                % print running speed
                app.SimulationResults.Value = [app.SimulationResults.Value; {'Simulation finished'}];
                app.SimulationResults.Value = [app.SimulationResults.Value; { ['Steady-state running speed is : ' num2str(Results.speed*3.6) 'km/h']}];
                app.SimulationResults.Value = [app.SimulationResults.Value; { 'Print score the results table'}];
                
                % plot results
                if isempty(ScoreVect)
                    ScoreVect = Results.speed*3.6;
                else
                    ScoreVect = [ScoreVect Results.speed*3.6];
                end
                plot(app.ScoreAxis,ScoreVect,'o','MarkerSize',6,'Color',[0 0 1],'MarkerFaceColor',[0 0 1]);
                drawnow; % to update all handles
            else % Print message to screen if simulation did not converge
                app.SimulationResults.Value = [app.SimulationResults.Value; {'Optimizer failed'}];
                app.SimulationResults.Value = [app.SimulationResults.Value; {'My worst nightmare just occured :), the optimization did not converge ...'}];
                app.SimulationResults.Value = [app.SimulationResults.Value; {'Please adjust the model parameters and execute the simulation again'}];
            end
            
            % print somethings to the screen
            %             SimulatonoutputTextAreaValueChanged()
        end

        % Value changed function: SimulationResults
        function SimulationResultsValueChanged(app, event)
            % figure ?
        end

        % Button pushed function: VideoButton
        function VideoButtonPushed(app, event)
            global S
            global vis
            global model_state;
            %global osimModel
            try
                motFile = S.MotFile;
                if ~exist(motFile,'file')
                    app.SimulationResults.Value = [app.SimulationResults.Value; {'Please run simulation first, using default simulation as visualisation'}];
                    motFile = fullfile(S.pathRepo,'Results',S.ResultsFolder,'DefaultSim_q.mot'); % file with outputs
                end
                dat = ReadMotFile(motFile);
                dt = dat.data(2,1)-dat.data(1,1);
                N = length(dat.data(:,1));
                
                CoordNamesAPI = {'pelvis_tilt','pelvis_tx','pelvis_ty','hip_flexion_r',...
                    'hip_flexion_l','lumbar_extension','knee_angle_r','knee_angle_l','ankle_angle_r','ankle_angle_l','mtp_angle_r','mtp_angle_l'};
                
                % column index in datastrcture for every coordName
                IndexCoord = nan(length(CoordNamesAPI),1);
                for i=1:length(CoordNamesAPI)
                    IndexCoord(i) = find(strcmp(CoordNamesAPI{i},dat.names))-1;
                end
                
                simbodyVis = vis.updSimbodyVisualizer();
                simbodyVis.setShowSimTime(true);
                Qvect = model_state.getQ();
                for i=[1:N 1]
                    dSel = dat.data(i,2:end);
                    for j=1:length(CoordNamesAPI)
                        if j==2 || j==3
                            Qvect.set(j-1,dSel(IndexCoord(j)));
                        else
                            Qvect.set(j-1,dSel(IndexCoord(j))*pi/180);
                        end
                    end
                    model_state.setQ(Qvect);
                    model_state.setTime(dat.data(i,1));
                    pause(dt);
                    vis.show(model_state);
                end
            catch
                app.SimulationResults.Value = [app.SimulationResults.Value; {'Unknown error in visualisation'}];
            end           
        end

        % Button pushed function: Button
        function ButtonPushed(app, event)
            global MForceAv
            if MForceAv>0
                MForceAv = MForceAv-50;
                app.MuscleForceNEditField.Value = MForceAv;
                app.SoleusEditField.Value = app.SoleusEditField.Value+50;
            else
                app.SimulationResults.Value = [app.SimulationResults.Value; {'You used all the additional muscle force'}];
            end
            %             drawnow
        end

        % Button pushed function: Button_2
        function Button_2Pushed(app, event)
            global MForceAv
            MForceAv = MForceAv+50;
            app.MuscleForceNEditField.Value = MForceAv;
            app.SoleusEditField.Value = app.SoleusEditField.Value-50;
            %             drawnow
        end

        % Button pushed function: Button_3
        function Button_3Pushed(app, event)
            global MForceAv
            if MForceAv>0
                MForceAv = MForceAv-50;
                app.MuscleForceNEditField.Value = MForceAv;
                app.GastrocnemiusEditField.Value = app.GastrocnemiusEditField.Value+50;
            else
                app.SimulationResults.Value = [app.SimulationResults.Value; {'You used all the additional muscle force'}];
            end
        end

        % Button pushed function: Button_4
        function Button_4Pushed(app, event)
            global MForceAv
            MForceAv = MForceAv+50;
            app.MuscleForceNEditField.Value = MForceAv;
            app.GastrocnemiusEditField.Value = app.GastrocnemiusEditField.Value-50;
        end

        % Button pushed function: Button_5
        function Button_5Pushed(app, event)
            global MForceAv
            if MForceAv>0
                MForceAv = MForceAv-50;
                app.MuscleForceNEditField.Value = MForceAv;
                app.RectusfemorisEditField.Value = app.RectusfemorisEditField.Value+50;
            else
                app.SimulationResults.Value = [app.SimulationResults.Value; {'You used all the additional muscle force'}];
            end
        end

        % Button pushed function: Button_6
        function Button_6Pushed(app, event)
            global MForceAv
            MForceAv = MForceAv+50;
            app.MuscleForceNEditField.Value = MForceAv;
            app.RectusfemorisEditField.Value = app.RectusfemorisEditField.Value-50;
        end

        % Button pushed function: ClearWindowButton
        function ClearWindowButtonPushed(app, event)
            app.SimulationResults.Value = {''};
        end

        % Value changed function: NormachillestendonstiffnessSlider
        function NormachillestendonstiffnessSliderValueChanged(app, event)
            kT = app.NormachillestendonstiffnessSlider.Value;
            
            % plot the tendon stiffness
            lTs = 0.2490;
            xV = 0.2490:0.001:0.27;
            lTtilde = xV./lTs;
            shift = getShift(kT);
            FV = (exp(kT.*(lTtilde - 0.995)))/5-0.25+shift;
            FV_default = (exp(35.*(lTtilde - 0.995)))/5-0.25+shift;
            FMo = 5137 + app.SoleusEditField.Value;
            plot(app.UIAxes,xV,FV_default.*FMo,'--k','LineWidth',1); hold(app.UIAxes, 'on');
            plot(app.UIAxes,xV,FV.*FMo,'b','LineWidth',2);
            hold(app.UIAxes, 'off')
            legend(app.UIAxes,{'Default','Adjusted'})
            drawnow
        end

        % Value changed function: FirstNameEditField
        function FirstNameEditFieldValueChanged(app, event)
            global S
            value = app.FirstNameEditField.Value;
            S.savename      = value;
        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            
            % add path to simulation code
            %TxtPath = which('MyID_PredSim2D_ab87ex.txt');
            %[MainPath,~,~] = fileparts(TxtPath);
            %rmpath(genpath(fullfile(MainPath)));
            
            % add the casadipath
            %CasFilePath = which('casadiMEX.mexw64');
            %[CasPath,~,~] = fileparts(CasFilePath);
            %rmpath(genpath(CasPath));
            
            % delete the app
            delete(app)
        end

        % Callback function
        function UIFigureCloseRequest2(app, event)
            delete(app)
            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1252 874];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create SimulateButton
            app.SimulateButton = uibutton(app.UIFigure, 'push');
            app.SimulateButton.ButtonPushedFcn = createCallbackFcn(app, @SimulateButtonPushed, true);
            app.SimulateButton.Position = [70 648 139 47];
            app.SimulateButton.Text = 'Simulate';

            % Create VideoButton
            app.VideoButton = uibutton(app.UIFigure, 'push');
            app.VideoButton.ButtonPushedFcn = createCallbackFcn(app, @VideoButtonPushed, true);
            app.VideoButton.Position = [458 648 150 42];
            app.VideoButton.Text = 'Video';

            % Create ResultsTableButton
            app.ResultsTableButton = uibutton(app.UIFigure, 'push');
            app.ResultsTableButton.Position = [259 650 150 42];
            app.ResultsTableButton.Text = 'Results Table';

            % Create SimulatonoutputTextAreaLabel
            app.SimulatonoutputTextAreaLabel = uilabel(app.UIFigure);
            app.SimulatonoutputTextAreaLabel.HorizontalAlignment = 'right';
            app.SimulatonoutputTextAreaLabel.FontSize = 14;
            app.SimulatonoutputTextAreaLabel.FontWeight = 'bold';
            app.SimulatonoutputTextAreaLabel.Position = [870 704 121 22];
            app.SimulatonoutputTextAreaLabel.Text = 'Simulaton output';

            % Create SimulationResults
            app.SimulationResults = uitextarea(app.UIFigure);
            app.SimulationResults.ValueChangedFcn = createCallbackFcn(app, @SimulationResultsValueChanged, true);
            app.SimulationResults.Position = [656 428 574 259];

            % Create FirstNameEditFieldLabel
            app.FirstNameEditFieldLabel = uilabel(app.UIFigure);
            app.FirstNameEditFieldLabel.HorizontalAlignment = 'right';
            app.FirstNameEditFieldLabel.Position = [70 754 64 22];
            app.FirstNameEditFieldLabel.Text = 'First Name';

            % Create FirstNameEditField
            app.FirstNameEditField = uieditfield(app.UIFigure, 'text');
            app.FirstNameEditField.ValueChangedFcn = createCallbackFcn(app, @FirstNameEditFieldValueChanged, true);
            app.FirstNameEditField.Position = [149 754 285 22];

            % Create DagvandewetenschapLabel
            app.DagvandewetenschapLabel = uilabel(app.UIFigure);
            app.DagvandewetenschapLabel.FontSize = 20;
            app.DagvandewetenschapLabel.FontWeight = 'bold';
            app.DagvandewetenschapLabel.Position = [270 822 396 42];
            app.DagvandewetenschapLabel.Text = 'Dag van de wetenschap: ';

            % Create Image
            app.Image = uiimage(app.UIFigure);
            app.Image.Position = [1080 764 100 100];
            app.Image.ImageSource = 'HMBLogo.png';

            % Create ClearWindowButton
            app.ClearWindowButton = uibutton(app.UIFigure, 'push');
            app.ClearWindowButton.ButtonPushedFcn = createCallbackFcn(app, @ClearWindowButtonPushed, true);
            app.ClearWindowButton.Position = [1110 708 100 22];
            app.ClearWindowButton.Text = 'Clear Window';

            % Create ScoreAxis
            app.ScoreAxis = uiaxes(app.UIFigure);
            title(app.ScoreAxis, 'Score History')
            xlabel(app.ScoreAxis, 'Attempt')
            ylabel(app.ScoreAxis, 'Running speed')
            app.ScoreAxis.Position = [37 428 554 160];

            % Create TendonstiffnessPanel
            app.TendonstiffnessPanel = uipanel(app.UIFigure);
            app.TendonstiffnessPanel.Title = 'Tendon stiffness';
            app.TendonstiffnessPanel.BackgroundColor = [0.7098 0.8549 0.9216];
            app.TendonstiffnessPanel.FontWeight = 'bold';
            app.TendonstiffnessPanel.FontSize = 14;
            app.TendonstiffnessPanel.Position = [656 15 584 378];

            % Create NormachillestendonstiffnessSliderLabel
            app.NormachillestendonstiffnessSliderLabel = uilabel(app.TendonstiffnessPanel);
            app.NormachillestendonstiffnessSliderLabel.HorizontalAlignment = 'right';
            app.NormachillestendonstiffnessSliderLabel.Position = [134 328 165 22];
            app.NormachillestendonstiffnessSliderLabel.Text = 'Norm achilles tendon stiffness';

            % Create NormachillestendonstiffnessSlider
            app.NormachillestendonstiffnessSlider = uislider(app.TendonstiffnessPanel);
            app.NormachillestendonstiffnessSlider.Limits = [20 50];
            app.NormachillestendonstiffnessSlider.ValueChangedFcn = createCallbackFcn(app, @NormachillestendonstiffnessSliderValueChanged, true);
            app.NormachillestendonstiffnessSlider.Position = [134 308 150 3];
            app.NormachillestendonstiffnessSlider.Value = 35;

            % Create UIAxes
            app.UIAxes = uiaxes(app.TendonstiffnessPanel);
            title(app.UIAxes, 'Achilles tendon stiffness')
            xlabel(app.UIAxes, 'Tendon length [m]')
            ylabel(app.UIAxes, 'Tendon Force [N]')
            app.UIAxes.BackgroundColor = [1 1 1];
            app.UIAxes.Position = [28 45 368 228];

            % Create Image3
            app.Image3 = uiimage(app.TendonstiffnessPanel);
            app.Image3.Position = [412 45 144 234];
            app.Image3.ImageSource = 'AchillesTendon.jpg';

            % Create MuscleforcePanel
            app.MuscleforcePanel = uipanel(app.UIFigure);
            app.MuscleforcePanel.Title = 'Muscle force ';
            app.MuscleforcePanel.BackgroundColor = [0.8588 0.8157 0.8157];
            app.MuscleforcePanel.FontWeight = 'bold';
            app.MuscleforcePanel.FontSize = 14;
            app.MuscleforcePanel.Position = [9 15 553 378];

            % Create Button
            app.Button = uibutton(app.MuscleforcePanel, 'push');
            app.Button.ButtonPushedFcn = createCallbackFcn(app, @ButtonPushed, true);
            app.Button.Position = [12 235 70 22];
            app.Button.Text = '+';

            % Create Button_2
            app.Button_2 = uibutton(app.MuscleforcePanel, 'push');
            app.Button_2.ButtonPushedFcn = createCallbackFcn(app, @Button_2Pushed, true);
            app.Button_2.Position = [12 205 70 22];
            app.Button_2.Text = '-';

            % Create Button_3
            app.Button_3 = uibutton(app.MuscleforcePanel, 'push');
            app.Button_3.ButtonPushedFcn = createCallbackFcn(app, @Button_3Pushed, true);
            app.Button_3.Position = [12 158 70 22];
            app.Button_3.Text = '+';

            % Create Button_4
            app.Button_4 = uibutton(app.MuscleforcePanel, 'push');
            app.Button_4.ButtonPushedFcn = createCallbackFcn(app, @Button_4Pushed, true);
            app.Button_4.Position = [12 128 70 22];
            app.Button_4.Text = '-';

            % Create Button_5
            app.Button_5 = uibutton(app.MuscleforcePanel, 'push');
            app.Button_5.ButtonPushedFcn = createCallbackFcn(app, @Button_5Pushed, true);
            app.Button_5.Position = [12 80 70 22];
            app.Button_5.Text = '+';

            % Create Button_6
            app.Button_6 = uibutton(app.MuscleforcePanel, 'push');
            app.Button_6.ButtonPushedFcn = createCallbackFcn(app, @Button_6Pushed, true);
            app.Button_6.Position = [12 50 70 22];
            app.Button_6.Text = '-';

            % Create MuscleForceNEditFieldLabel
            app.MuscleForceNEditFieldLabel = uilabel(app.MuscleforcePanel);
            app.MuscleForceNEditFieldLabel.HorizontalAlignment = 'right';
            app.MuscleForceNEditFieldLabel.FontSize = 14;
            app.MuscleForceNEditFieldLabel.Position = [-3 302 111 22];
            app.MuscleForceNEditFieldLabel.Text = 'Muscle Force [N]';

            % Create MuscleForceNEditField
            app.MuscleForceNEditField = uieditfield(app.MuscleforcePanel, 'numeric');
            app.MuscleForceNEditField.Position = [123 302 100 22];
            app.MuscleForceNEditField.Value = 1000;

            % Create SoleusEditFieldLabel
            app.SoleusEditFieldLabel = uilabel(app.MuscleforcePanel);
            app.SoleusEditFieldLabel.HorizontalAlignment = 'right';
            app.SoleusEditFieldLabel.FontSize = 14;
            app.SoleusEditFieldLabel.FontWeight = 'bold';
            app.SoleusEditFieldLabel.Position = [140 235 51 22];
            app.SoleusEditFieldLabel.Text = 'Soleus';

            % Create SoleusEditField
            app.SoleusEditField = uieditfield(app.MuscleforcePanel, 'numeric');
            app.SoleusEditField.Position = [120 211 100 22];

            % Create GastrocnemiusEditFieldLabel
            app.GastrocnemiusEditFieldLabel = uilabel(app.MuscleforcePanel);
            app.GastrocnemiusEditFieldLabel.HorizontalAlignment = 'right';
            app.GastrocnemiusEditFieldLabel.FontSize = 14;
            app.GastrocnemiusEditFieldLabel.FontWeight = 'bold';
            app.GastrocnemiusEditFieldLabel.Position = [114 158 107 22];
            app.GastrocnemiusEditFieldLabel.Text = 'Gastrocnemius';

            % Create GastrocnemiusEditField
            app.GastrocnemiusEditField = uieditfield(app.MuscleforcePanel, 'numeric');
            app.GastrocnemiusEditField.Position = [123 137 100 22];

            % Create RectusfemorisEditFieldLabel
            app.RectusfemorisEditFieldLabel = uilabel(app.MuscleforcePanel);
            app.RectusfemorisEditFieldLabel.HorizontalAlignment = 'right';
            app.RectusfemorisEditFieldLabel.FontSize = 14;
            app.RectusfemorisEditFieldLabel.FontWeight = 'bold';
            app.RectusfemorisEditFieldLabel.Position = [112 80 107 22];
            app.RectusfemorisEditFieldLabel.Text = 'Rectus femoris';

            % Create RectusfemorisEditField
            app.RectusfemorisEditField = uieditfield(app.MuscleforcePanel, 'numeric');
            app.RectusfemorisEditField.Position = [118 59 100 22];

            % Create MuscleForceLabel
            app.MuscleForceLabel = uilabel(app.MuscleforcePanel);
            app.MuscleForceLabel.FontSize = 14;
            app.MuscleForceLabel.FontWeight = 'bold';
            app.MuscleForceLabel.Position = [123 324 144 36];
            app.MuscleForceLabel.Text = 'Muscle Force';

            % Create Image2
            app.Image2 = uiimage(app.MuscleforcePanel);
            app.Image2.Position = [276 5 243 368];
            app.Image2.ImageSource = '2DModel_LabelMuscles-01.png';

            % Create ReadyLampLabel
            app.ReadyLampLabel = uilabel(app.UIFigure);
            app.ReadyLampLabel.HorizontalAlignment = 'right';
            app.ReadyLampLabel.Position = [673 704 50 22];
            app.ReadyLampLabel.Text = 'Ready ?';

            % Create ReadyLamp
            app.ReadyLamp = uilamp(app.UIFigure);
            app.ReadyLamp.Position = [726 700 30 30];

            % Create MusculoskeletalmodelLabel
            app.MusculoskeletalmodelLabel = uilabel(app.UIFigure);
            app.MusculoskeletalmodelLabel.FontSize = 16;
            app.MusculoskeletalmodelLabel.FontWeight = 'bold';
            app.MusculoskeletalmodelLabel.Position = [825 804 191 42];
            app.MusculoskeletalmodelLabel.Text = 'Musculoskeletal model';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = OptSpeedRunner_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end