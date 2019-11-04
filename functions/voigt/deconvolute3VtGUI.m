function deconvolute3VtGUI()

global x y xInit yInit
% deconvolute3VtGUI   Deconvolution of the Mordenite zeolite spectrum by 3 Voigt bands 


% Main window
hwin = figure('Visible','off','Position',[50 50 900 500],'NumberTitle','off','Name','Mordenite IR spectra deconvolution',...
    'MenuBar','none','Toolbar','none','Resize','on','ResizeFcn',{@WinResizeFcn});
bgclr = get(hwin,'Color');

% Axes
haxes = axes('Units','pixels','XDir','reverse');
xlabel(haxes,'Wavenumber, cm^-^1');
ylabel(haxes,'Absorbance, a.u.');
% Components
hLoad = uicontrol('Style','pushbutton','String','Load','Position',[20,400,140,25],'Callback',{@LoadFcn},'BackgroundColor',bgclr); 
hBackground = uicontrol('Style','pushbutton','String','Remove Background','Position',[20,360,140,25],'Callback',{@BackFcn},'BackgroundColor',bgclr,'Enable','off'); 
hDeconvolute = uicontrol('Style','pushbutton','String','Deconvolution','Position',[20,320,140,25],'Callback',{@DeconvFcn},'BackgroundColor',bgclr,'Enable','off'); 

hLoadText = uicontrol('Style','text','String','Load spectrum from a file. Only 2-column XY file format is supported','HorizontalAlignment','left','BackgroundColor',bgclr);
hCopyText = uicontrol('Style','text','String','See DOI:XXXXXXXXX for details','HorizontalAlignment','left','BackgroundColor',bgclr);
hText = uicontrol('Style','text','String','Deconvolutions:','HorizontalAlignment','left','BackgroundColor',bgclr);
hDeconvText = uicontrol('Style','edit','String','1','BackgroundColor','white');
hOutText = uicontrol('Style','text','String','','HorizontalAlignment','left','BackgroundColor',bgclr);

% Move the GUI to the center of the screen
movegui(hwin,'center');

% Make the GUI visible
set(hwin,'Visible','on');

% Callback functions
    function WinResizeFcn(source,eventdata)
        % Resize the window
        pos = get(hwin,'Position');
        border=pos(1:2);
        w = pos(3);
        h = pos(4);
        if w<600 || h<450, set(hwin,'Position',[border 600 450]);end
        
        set(hLoadText,'Position',   [20,    h-30,140,25]);
        set(hLoad,'Position',       [20,    h-60,140,25]);
        set(hBackground,'Position', [20,    h-90,140,25]);
        set(hText,'Position',       [20,    h-120,140,25]);
        set(hDeconvText,'Position', [20,    h-140,140,25]);
        set(hDeconvolute,'Position',[20,    h-170,140,25]);
        set(hOutText,'Position',    [20,    h-360,160,155]);
        
        set(hCopyText,'Position',   [10,    10,160,17]);
        set(haxes,'Position',       [240,   50   ,w-260,h-70]);
    end

    function LoadFcn(source,eventdata)
        % Load the spectrum in XY format
        [filename, pathname] = uigetfile('*.txt;*.csv;*.xy','Load the spectrum file in XY file format');
        
        if (filename~=0),
            loaddata = fullfile(pathname,filename);
            
            % the file must contain only TWO columns, X and Y - nothing else
            data = load(loaddata,'-ascii');
            x=data(:,1);
            y=data(:,2);
        end;
        plot(haxes,x,y);
        set(haxes,'XDir','reverse');
        xlabel(haxes,'Wavenumber, cm^-^1');
        ylabel(haxes,'Absorbance, a.u.');
        xInit=x;yInit=y;
        set(hBackground,'Enable','on');
        set(hDeconvolute,'Enable','off');
        set(hOutText,'String','');
    end
  
    function BackFcn(source,eventdata)
        % Removes background 
        sel=(x>3500 & x<3700);
        x=x(sel);y=y(sel);
        
        % indexes of the points closest to 3500 and 3700 1/cm
        [~, i1]=min(abs(x-3500));[~, i2]=min(abs(x-3700));
        p=polyfit([x(i1) x(i2)],[y(i1) y(i2)],1);
        y=y-polyval(p,x);
        plot(haxes,x,y);
        set(haxes,'XDir','reverse');
        xlabel(haxes,'Wavenumber, cm^-^1');
        ylabel(haxes,'Absorbance, a.u.');
        yLimit = get(haxes,'ylim');ylim([0 yLimit(2)]);
        set(hDeconvolute,'Enable','on');
    end
    
    function DeconvFcn(source,eventdata)
        % deconvolution of the spectrum
        
        % if you have read so and have the Parralel Computation Toolbox,
        % you should run the following command MATLAB 
        % matlabpool('open',XXX),
        % here XXX - number of real cores in your CPU
        nIterations = get(hDeconvText,'String');nIterations = str2double(nIterations);
        [fittingError areaBands plots]=deconvolute3Vt(x,y,nIterations);
        
        % drawing the deconvolution
        plot(haxes,x,y,'-k');
        hold on;                
        colors={'g','r','b','m'};
        for ii=1:4,
            plot(haxes,x,plots(:,ii),colors{ii});
            plot(haxes,x,plots(:,ii+4),colors{ii});
            fill([x' fliplr(x')],[plots(:,ii)' fliplr(plots(:,ii+4)')],colors{ii},'EdgeColor','none');alpha(.4);
        end;
        hold off;
        yLimit = get(haxes,'ylim');ylim([0 yLimit(2)]);
        set(haxes,'XDir','reverse');
        xlabel(haxes,'Wavenumber, cm^-^1');
        ylabel(haxes,'Absorbance, a.u.');
        
        outString=sprintf('Fitting error: %4.1f \t-\t%4.1f\t%%\n\nArea fractions:\nLF = %4.1f\t-\t%4.1f\t%%\nHF = %4.1f\t-\t%4.1f \t%%\nTF = %4.1f\t-\t%4.1f \t%%\n\n',fittingError(1),fittingError(2),areaBands(1),areaBands(2),areaBands(3),areaBands(4),areaBands(5),areaBands(6));
        extinction=[1 0.73 0.62];   % LF, HF, and TF relative extinction coefficients
        BASfraction=zeros(6,1);
        for ii=1:3,
            BASfraction(ii)=    areaBands(ii*2-1)   /extinction(ii)/(areaBands(2)/extinction(1)+areaBands(4)/extinction(2)+areaBands(6)/extinction(3))*100;      %min
            BASfraction(ii+3)=  areaBands(ii*2)	/extinction(ii)/(areaBands(1)/extinction(1)+areaBands(3)/extinction(2)+areaBands(5)/extinction(3))*100;      %max
        end;
        outString=[outString,sprintf('BAS fractions: \nLF = %4.1f\t-\t%4.1f\t%%\nHF = %4.1f\t-\t%4.1f \t%%\nTF = %4.1f\t-\t%4.1f \t%%',BASfraction(1),BASfraction(4),BASfraction(2),BASfraction(5),BASfraction(3),BASfraction(6))];
        set(hOutText,'String',outString);
    end


    



uiwait(gcf);

end