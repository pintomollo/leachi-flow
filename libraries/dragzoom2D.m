function dragzoom2D(hAx)
%DRAGZOOM2D Drag and zoom tool simplified from dragzoom
%
% Description:
%   DRAGZOOM2D allows you to interactively manage the axes in figure.
%   This simple function for usable draging and zooming of axes, using the 
%   mouse and the keyboard shortcuts.
%
% Using:
%   dragzoom(hAx)
%
% Interactive mode:
%   Mouse actions:
%       single-click and holding LB : Activation Drag mode
%       single-click and holding RB : Activation Rubber Band for region zooming
%       single-click MB             : Activation 'Extend' Zoom mode
%       double-click LB, RB, MB     : Reset to Original View
% 
%   Hotkeys:
%       '+'                         : Zoom plus
%       '-'                         : Zoom minus
%       '0'                         : Set default axes (reset to original view)
%       'uparrow'                   : Up or down (inrerse) drag
%       'downarrow'                 : Down or up (inverse) drag
%       'leftarrow'                 : Left or right (inverse) drag
%       'rightarrow'                : Right or left (inverse) drag
%       'x'                         : If pressed, toggles zoom and drag only for X axis
%       'y'                         : If pressed, toggles zoom and drag only for Y axis
%

    if (nargin == 1 && ishandle(hAx) && strncmp(get(hAx, 'Type'), 'axes', 4))
        hFig = ancestor(hAx, 'figure');
    else
        error('dragzoom2D:invalidInputs', ...
            'Input must be an axes handle.')
    end

    params = Setup(hAx);

    set(hFig, 'CurrentAxes', hAx, ...
        'userdata', params, ...
        'WindowButtonDownFcn',      {@WindowButtonDownCallback2D}, ...
        'WindowButtonUpFcn',        {@WindowButtonUpCallback2D}, ...
        'WindowButtonMotionFcn',    {@WindowButtonMotionCallback2D}, ...
        'KeyPressFcn',              {@WindowKeyPressCallback2D}, ...
        'KeyReleaseFcn',            {@WindowKeyReleaseCallback2D}, ...
        'WindowKeyPressFcn',        {@WindowKeyPressCallback2D}, ...
        'WindowKeyReleaseFcn',      {@WindowKeyReleaseCallback2D});

    return;
end

%--------------------------------------------------------------------------

%==========================================================================
function WindowButtonDownCallback2D(src, evnt)    %#ok
    %WindowButtonDownCallback2D
    
    params = get(src, 'userdata');
    clickType = get(src, 'SelectionType');
    
    switch clickType
        case 'normal'
            %----DragMouseBegin();
            if (~params.fIsDragAllowed)
                [cx, cy] = GetCursorCoordOnWindow(src);
                
                params.mStartX = cx;
                params.mStartY = cy;
                
                params.fIsDragAllowed = true;
            end
        case 'open'
            %---ResetAxesToOrigView();
            SetAxesLimits(params.hAx, params.mDefaultXLim, params.mDefaultYLim);
            params.mZoomIndexX = find(params.mZoomGrid == 100);
            params.mZoomIndexY = params.mZoomIndexX;
        case 'alt'
            %---RubberBandBegin();
            if (~params.fIsRubberBandOn)
                [acx, acy] = GetCursorCoordOnAxes(params.hAx);
               
                %---RubberBandSetup();
                h1 = patch([acx acx acx acx], [acy acy acy acy], 'k',
                    'Parent', params.hAx, ...
                    'EdgeColor', 'w', ...
                    'FaceColor', 'none', ...
                    'LineWidth', 1.5, ...
                    'LineStyle', '-');
                
                h2 = patch([acx acx acx acx], [acy acy acy acy], 'k',
                    'Parent', params.hAx, ...
                    'EdgeColor', params.mRbEdgeColor, ...
                    'FaceColor', params.mRbFaceColor, ...
                    'FaceAlpha', params.mRbFaceAlpha, ...
                    'LineWidth', 0.5, ...
                    'LineStyle', '-');    
 
                % create rubber band struct
                params.mRubberBand = struct(...
                    'obj',	[h1 h2], ...
                    'x1',  	acx, ...
                    'y1',  	acy, ...
                    'x2',  	acx, ...
                    'y2',  	acy);
                
                params.fIsRubberBandOn = true;
            end
        case 'extend'
            %---ZoomMouseExtendBegin();
            if ~params.fIsZoomExtendAllowed
                %UpdateCurrentZoomAxes();
                
                % set new zoom grid for extend zoom
                [params.mZoomGrid, params.mZoomSteps] = ZoomLogGrid(params.mZoomMinPow, params.mZoomMaxPow, params.mZoomExtendNum);
                
                %---UpdateCurrentZoomAxes();
                [xLim, yLim] = GetAxesLimits(params.hAx);
                [curentZoomX, curentZoomY] = GetCurrentZoomAxesPercent(params.hAx, xLim, yLim, params.mDefaultXLim, params.mDefaultYLim);
                [nu, params.mZoomIndexX] = min(abs(params.mZoomGrid - curentZoomX));  %#ok ([~, ...])
                [nu, params.mZoomIndexY] = min(abs(params.mZoomGrid - curentZoomY));  %#ok ([~, ...])

                [wcx, wcy] = GetCursorCoordOnWindow(src);
                [acx, acy] = GetCursorCoordOnAxes(params.hAx);
                
                params.mStartX = wcx;
                params.mStartY = wcy;

                params.mBindX = acx;
                params.mBindY = acy;

                params.fIsZoomExtendAllowed = true;
            end
    end

    set(src, 'userdata', params);
end
%--------------------------------------------------------------------------

%==========================================================================
function WindowButtonUpCallback2D(src, evnt)      %#ok
    %WindowButtonUpCallback2D
    
    params = get(src, 'userdata');

    %---DragMouseEnd();
    if params.fIsDragAllowed
        params.fIsDragAllowed = false;
    end

    %---ZoomMouseExtendEnd();
    if params.fIsZoomExtendAllowed
        % set default zoom grid
        %---SetDefaultZoomGrid();
        params.mZoomGrid = params.mDefaultZoomGrid;
        params.mZoomSteps = params.mDefaultZoomSteps;
        
        params.mZoomIndexX = params.find(mZoomGrid == 100);
        params.mZoomIndexY = params.mZoomIndexX;

        params.fIsZoomExtendAllowed = false;
    end

    %---RubberBandEnd();
    if params.fIsRubberBandOn
        params.fIsRubberBandOn = false;
        
        delete(params.mRubberBand.obj);          
        %---RubberBandZoomAxes();
        xLim = sort([params.mRubberBand.x1, params.mRubberBand.x2]);
        yLim = sort([params.mRubberBand.y1, params.mRubberBand.y2]);
        
        if (range(xLim) ~= 0 && range(yLim) ~= 0)
            [zoomPctX, zoomPctY] = GetCurrentZoomAxesPercent(params.hAx, ...
                xLim, yLim, params.mDefaultXLim, params.mDefaultYLim);
            
            if params.fIsImage
                zoomPctX = min(zoomPctX, zoomPctY);
                zoomPctY = zoomPctX;
            end
            
            cx = mean(xLim);
            cy = mean(yLim);
            
            xLim = RecalcZoomAxesLimits(xLim, params.mDefaultXLim, ...
                cx, zoomPctX, strcmp(get(params.hAx, 'xscale'), 'log'));
            yLim = RecalcZoomAxesLimits(yLim, params.mDefaultYLim, ...
                cy, zoomPctY, strcmp(get(params.hAx, 'yscale'), 'log'));
            
            SetAxesLimits(params.hAx, xLim, yLim);
        end

        params.mRubberBand = [];
    end

    set(src, 'userdata', params);
end
%--------------------------------------------------------------------------

%==========================================================================
function WindowButtonMotionCallback2D(src, evnt)  %#ok
    %WindowButtonMotionCallback2D
            
    params = get(src, 'userdata');

    %---DragMouse();
    if params.fIsDragAllowed
        [cx, cy] = GetCursorCoordOnWindow(src);
        
        pdx = params.mStartX - cx;
        pdy = params.mStartY - cy;
        
        params.mStartX = cx;
        params.mStartY = cy;

        if (params.fIsImage)
            pdy = -pdy;
        end
        
        DragAxes(params.hAx, pdx, pdy, ...
            params.fIsEnableDragX, params.fIsEnableDragY);
    end

    %---RubberBandUpdate();
    if params.fIsRubberBandOn
        [acx, acy] = GetCursorCoordOnAxes(params.hAx);
        
        params.mRubberBand.x2 = acx;
        params.mRubberBand.y2 = acy;

        %---RubberBandSetPos();
        x1 = params.mRubberBand.x1;
        y1 = params.mRubberBand.y1;
        
        set(params.mRubberBand.obj, ...
            'XData', [x1 acx acx x1], ...
            'YData', [y1 y1 acy acy]);
    end
    
    %---ZoomMouseExtend();
    if params.fIsZoomExtendAllowed
        
        % Heuristic for pixel change to camera zoom factor 
        % (taken from function ZOOM)
        [wcx, wcy] = GetCursorCoordOnWindow(src);
        
        xy(1) = wcx - params.mStartX;
        xy(2) = wcy - params.mStartY;
        q = max(-0.9, min(0.9, sum(xy)/70)) + 1;

        %directions = {'minus', 'plus'};
        if (q < 1)
            %direction = directions{1};
            dz = -1;
        elseif (q > 1)
            %direction = directions{2};
            dz = -1;
        else
            return;
        end
        
        %---ZoomAxes(direction, mZoom3DBindX, mZoom3DBindY)
        [xLim, yLim] = GetAxesLimits(params.hAx);
        
        if params.fIsImage
            params.mZoomIndexX = params.mZoomIndexX + dz;
            params.mZoomIndexY = params.mZoomIndexY;
        else
            if params.fIsEnableZoomX
                params.mZoomIndexX = params.mZoomIndexX + dz;
            end
            if params.fIsEnableZoomY
                params.mZoomIndexY = params.mZoomIndexY + dz;
            end
        end
        xLim = RecalcZoomAxesLimits(xLim, params.mDefaultXLim, ...
            params.mBindX, params.mZoomGrid(params.mZoomIndexX), ...
            strcmp(get(params.hAx, 'xscale'), 'log'));
        yLim = RecalcZoomAxesLimits(yLim, params.mDefaultYLim, ...
            params.mBindY, params.mZoomGrid(params.mZoomIndexY), ...
            strcmp(get(params.hAx, 'yscale'), 'log'));
        
        SetAxesLimits(params.hAx, xLim, yLim);
        
        params.mStartX = wcx;
        params.mStartY = wcy;            
    end

    set(src, 'userdata', params);
end
%--------------------------------------------------------------------------

%==========================================================================
function WindowKeyPressCallback2D(src, evnt)      %#ok
    %WindowKeyPressCallback2D

    params = get(src, 'userdata');
    modifier = evnt.Modifier;
    dz = 0;

    switch evnt.Key
        case '0'
            %---ResetAxesToOrigView();
            SetAxesLimits(params.hAx, params.mDefaultXLim, params.mDefaultYLim);
            params.mZoomIndexX = find(params.mZoomGrid == 100);
            params.mZoomIndexY = params.mZoomIndexX;
        case {'equal', 'add', 'plus', '1'}
            %ZoomKeys('plus');
            dz = 1;
        case {'hyphen', 'subtract', 'minus', 'slash'}
            %ZoomKeys('minus');
            dz = -1;
        case {'leftarrow', 'rightarrow', 'uparrow', 'downarrow', ...}
                'left', 'right', 'up', 'down'}

            %---DragKeys('...');
            dx = params.mDragShiftStep;
            dy = params.mDragShiftStep;
            
            % Increment of speed when you hold the button
            params.mDragShiftStep = params.mDragShiftStep + params.mDragShiftStepInc;
            
            switch evnt.Key(1)
                case 'r'
                    dx = -dx;
                    dy = 0;
                case 'l'
                    dy = 0;
                case 'd'
                    dx = 0;
                case 'u'
                    dx = 0;
                    dy = -dy;
            end

            if (params.fIsImage)
                pdy = -pdy;
            end

            DragAxes(params.hAx, dx, dy, ...
                params.fIsEnableDragX, params.fIsEnableDragY);
    
        case 'x'
            %DragEnable('y', 'off');
            %ZoomEnable('y', 'off');
            params.fIsEnableDragY = ~params.fIsEnableDragY;
            params.fIsEnableZoomY = ~params.fIsEnableZoomY;
        case 'y'
            %DragEnable('x', 'off');
            %ZoomEnable('x', 'off');
            params.fIsEnableDragX = ~params.fIsEnableDragX;
            params.fIsEnableZoomX = ~params.fIsEnableZoomX;
    end

    if (dz ~= 0)
        %UpdateCurrentZoomAxes();
        
        % set new zoom grid for extend zoom
        [params.mZoomGrid, params.mZoomSteps] = ZoomLogGrid(params.mZoomMinPow, params.mZoomMaxPow, params.mZoomExtendNum);
        
        %---UpdateCurrentZoomAxes();
        [xLim, yLim] = GetAxesLimits(params.hAx);
        [curentZoomX, curentZoomY] = GetCurrentZoomAxesPercent(params.hAx, xLim, yLim, params.mDefaultXLim, params.mDefaultYLim);
        [nu, params.mZoomIndexX] = min(abs(params.mZoomGrid - curentZoomX));  %#ok ([~, ...])
        [nu, params.mZoomIndexY] = min(abs(params.mZoomGrid - curentZoomY));  %#ok ([~, ...])

        [acx, acy] = GetCursorCoordOnAxes(params.hAx);
        
        if params.fIsImage
            params.mZoomIndexX = params.mZoomIndexX + dz;
            params.mZoomIndexY = params.mZoomIndexX;
        else
            if params.fIsEnableZoomX
                params.mZoomIndexX = params.mZoomIndexX + dz;
            end
            if params.fIsEnableZoomY
                params.mZoomIndexY = params.mZoomIndexY + dz;
            end
        end
        xLim = RecalcZoomAxesLimits(xLim, params.mDefaultXLim, ...
            acx, params.mZoomGrid(params.mZoomIndexX), ...
            strcmp(get(params.hAx, 'xscale'), 'log'));
        yLim = RecalcZoomAxesLimits(yLim, params.mDefaultYLim, ...
            acy, params.mZoomGrid(params.mZoomIndexY), ...
            strcmp(get(params.hAx, 'yscale'), 'log'));
        
        SetAxesLimits(params.hAx, xLim, yLim);
        
        %ZoomAxes(direction, acx, acy)
        %PointerCrossUpdate();
        %SetDefaultZoomGrid();
    end

    set(src, 'userdata', params);
end
%--------------------------------------------------------------------------

%==========================================================================
function WindowKeyReleaseCallback2D(src, evnt)    %#ok
    %WindowKeyReleaseCallback2D
    
    params = get(src, 'userdata');

    switch evnt.Key
        case {'leftarrow', 'rightarrow', 'uparrow', 'downarrow'}
            params.mDragShiftStep = params.mDragSaveShiftStep;
        %case 'x'
            %DragEnable('y', 'on');
            %ZoomEnable('y', 'on');
        %    params.fIsEnableDragY = true;
        %    params.fIsEnableZoomY = true;
        %case 'y'
            %DragEnable('x', 'on');
            %ZoomEnable('x', 'on');
        %    params.fIsEnableDragY = true;
        %    params.fIsEnableZoomY = true;
    end
end
%--------------------------------------------------------------------------

%==========================================================================
function params = Setup(hAx)
    %Setup setup options
            
    h = findobj(hAx, 'Type', 'Image');
        
    mZoomMinPow = 0;
    mZoomMaxPow = 5;
    mZoomNum = 51;
    mZoomIndex = 11; % 100%
    [mDefaultZoomGrid, mDefaultZoomSteps] = ...
        ZoomLogGrid(mZoomMinPow, mZoomMaxPow, mZoomNum);
    
    % handles
    params = struct('hAx', hAx, ...
        'mStartX', [], ...
        'mStartY', [], ...
        'mBindX', [], ...
        'mBindY', [], ...
        'mDragShiftStep', 3, ...
        'mDragSaveShiftStep', 3, ...
        'mDragShiftStepInc', 1, ...
        'mZoomMinPow', mZoomMinPow, ...
        'mZoomMaxPow', mZoomMaxPow, ...
        'mZoomNum', mZoomNum, ...
        'mZoomExtendNum', 301, ...
        'mZoomKeysNum', 181, ...
        'mDefaultZoomGrid', mDefaultZoomGrid, ...
        'mDefaultZoomSteps', mDefaultZoomSteps, ...
        'mZoomGrid', mDefaultZoomGrid, ...
        'mZoomSteps', mDefaultZoomSteps, ...
        'mZoomIndexX', mZoomIndex, ...
        'mZoomIndexY', mZoomIndex, ...
        'mDefaultXLim', get(hAx, 'XLim'), ...
        'mDefaultYLim', get(hAx, 'YLim'), ...
        'mRubberBand', [], ...
        'mRbEdgeColor', 'k', ...
        'mRbFaceColor', 'none', ...
        'mRbFaceAlpha', 1, ...
        'fIsDragAllowed', false, ...
        'fIsZoomExtendAllowed', false, ...
        'fIsRubberBandOn', false, ...
        'fIsEnableDragX', true, ...
        'fIsEnableDragY', true, ...
        'fIsEnableZoomX', true, ...
        'fIsEnableZoomY', true, ...
        'fIsImage', ~isempty(h));

end
%--------------------------------------------------------------------------

%==========================================================================
function DragAxes(hAx, pdx, pdy, fIsEnableDragX, fIsEnableDragY)
    %DragAxes
    
    [xLim, yLim] = GetAxesLimits(hAx);
    
    pos = GetObjPos(hAx, 'Pixels');
    pbar = get(hAx, 'PlotBoxAspectRatio');
    
    %NOTE: MATLAB Bug?
    % Fixed problem with AspectRatio and Position of Axes
    % MATLAB Function PAN is not correct works with rectangular images!
    % Here it is correctly.
    
    imAspectRatioX = pbar(2) / pbar(1);
    if (imAspectRatioX ~= 1)
        posAspectRatioX = pos(3) / pos(4);
        arFactorX = imAspectRatioX * posAspectRatioX;
        if (arFactorX < 1)
            arFactorX = 1;
        end
    else
        arFactorX = 1;
    end
    
    imAspectRatioY = pbar(1) / pbar(2);
    if (imAspectRatioY ~= 1)
        posAspectRatioY = pos(4) / pos(3);
        arFactorY = imAspectRatioY * posAspectRatioY;
        if (arFactorY < 1)
            arFactorY = 1;
        end
    else
        arFactorY = 1;
    end
    
    if fIsEnableDragX
        % For log plots, transform to linear scale
        if strcmp(get(hAx, 'xscale'), 'log')
            xLim = log10(xLim);
            xLim = FixInfLogLimits('x', xLim);
            isXLog = true;
        else
            isXLog = false;
        end
        
        dx = pdx * range(xLim) / (pos(3) / arFactorX);
        xLim = xLim + dx;
        
        % For log plots, untransform limits
        if isXLog
            xLim = 10.^(xLim);
        end
    end
    if fIsEnableDragY
        if strcmp(get(hAx, 'yscale'), 'log')
            yLim = log10(yLim);
            yLim = FixInfLogLimits('y', yLim);
            isYLog = true;
        else
            isYLog = false;
        end
        
        dy = pdy * range(yLim) / (pos(4) / arFactorY);
        
        yLim = yLim + dy; 
        
        if isYLog
            yLim = 10.^(yLim);
        end
    end
    
    SetAxesLimits(hAx, xLim, yLim);
end
%--------------------------------------------------------------------------

%==========================================================================
function [curentZoomX, curentZoomY] = GetCurrentZoomAxesPercent(hAx, xLim, yLim, mDefaultXLim, mDefaultYLim)
    %GetCurrentZoomAxesPercent
    
    if strcmp(get(hAx, 'xscale'), 'log')
        xLim = log10(xLim);
        defaultXLim = log10(mDefaultXLim);
    else
        defaultXLim = mDefaultXLim;
    end
    if strcmp(get(hAx, 'yscale'), 'log')
        yLim = log10(yLim);
        defaultYLim = log10(mDefaultYLim);
    else
        defaultYLim = mDefaultYLim;
    end
    
    curentZoomX = range(defaultXLim) * 100 / range(xLim);
    curentZoomY = range(defaultYLim) * 100 / range(yLim);
end
%--------------------------------------------------------------------------

%==========================================================================
    function [x, y, z] = GetCursorCoordOnAxes(hAx)
        %GetCursorCoordOnAxImg
        
        crd = get(hAx, 'CurrentPoint');
        x = crd(2,1);
        y = crd(2,2);
        z = crd(2,3);
    end
%--------------------------------------------------------------------------

%==========================================================================
    function [x, y] = GetCursorCoordOnWindow(hFig)
        %GetCursorCoordOnWindow
        
        dfltUnits = get(hFig, 'Units');
        set(hFig, 'Units', 'pixels');
        
        crd = get(hFig, 'CurrentPoint');
        x = crd(1); 
        y = crd(2);
        
        set(hFig, 'Units', dfltUnits);
    end
%--------------------------------------------------------------------------

%==========================================================================
function pos = GetObjPos(h, units)
    %GetObjPos get object position
    
    dfltUnits = get(h, 'Units');
    set(h, 'Units', units);
    pos = get(h, 'Position');
    set(h, 'Units', dfltUnits);
end
%--------------------------------------------------------------------------

%==========================================================================
    function [xLim, yLim] = GetAxesLimits(hAx)
        %GetAxesLimits
        
        xLim = get(hAx, 'XLim');
        yLim = get(hAx, 'YLim');
    end
%--------------------------------------------------------------------------

%==========================================================================
    function SetAxesLimits(hAx, xLim, yLim)
        %SetAxesLimits
        
        set(hAx, 'XLim', xLim);
        set(hAx, 'YLim', yLim);
    end
%--------------------------------------------------------------------------

%==========================================================================
    function axLim = RecalcZoomAxesLimits(axLim, axLimDflt, zcCrd, zoomPct, isLog)
        %RecalcZoomAxesLimits recalc axes limits
        
        if isLog
            axLim = log10(axLim);
            %axLim = FixInfLogLimits(ax, axLim);
            axLimDflt = log10(axLimDflt);
            zcCrd = log10(zcCrd);

            if (~all(isfinite(axLim)) || ~all(isreal(axLim)))
                axLim = axLimDflt;
            end
        end
                
        if (zcCrd < axLim(1)), zcCrd = axLim(1); end
        if (zcCrd > axLim(2)), zcCrd = axLim(2); end
        
        rf = range(axLim);
        ra = range([axLim(1), zcCrd]);
        rb = range([zcCrd, axLim(2)]);
        
        cfa = ra / rf; 
        cfb = rb / rf;
        
        newRange = range(axLimDflt) * 100 / zoomPct;
        dRange = newRange - rf;
        
        axLim(1) = axLim(1) - dRange * cfa;
        axLim(2) = axLim(2) + dRange * cfb;
        
        if isLog
            axLim = 10.^axLim;
        end
    end
%--------------------------------------------------------------------------

%==========================================================================
    function [zg, st] = ZoomLogGrid(a, b, n)
        %ZoomLogGrid log zoom grid
        
        zg = unique(round(logspace(a, b, n)));
        
        zg(zg<10) = [];	% begin zoom == 10%
        st = length(zg);
        
    end
%--------------------------------------------------------------------------

