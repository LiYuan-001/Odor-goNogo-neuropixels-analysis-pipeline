function eventTable = annotateSwitchEvents(switchSig, Fs)
% annotateSwitchEvents_anyClick
% Interactive annotation of switch signal events
%
% Mouse:
%   Left click   -> assign event at clicked time
%   Right click  -> set x-axis to [x-20, x+20]
%   Scroll       -> zoom x-axis in/out around current mouse x-position
%
% Buttons:
%   Finish     -> stop interaction, print table, return from function,
%                 keep figure open for viewing
%   Clear last -> remove last assigned event
%   Clear all  -> remove all assigned events
%   Reset X    -> restore original x-axis range
%   Close      -> close figure
%
% Output:
%   eventTable with columns:
%   - eventName
%   - signalIdx

if nargin < 2 || isempty(Fs)
    Fs = 1;
end

switchSig = switchSig(:);
nSamp = length(switchSig);
t = (0:nSamp-1)' / Fs;

% Make binary if needed
u = unique(switchSig(~isnan(switchSig)));
if ~all(ismember(u, [0 1]))
    th = (max(switchSig) + min(switchSig)) / 2;
    switchSig = switchSig > th;
end

maxEvt = 10;

% State
S.switchSig = switchSig;
S.Fs = Fs;
S.nSamp = nSamp;
S.t = t;
S.maxEvt = maxEvt;

S.activeBox = 1;
S.selectedIdx = nan(maxEvt,1);

S.editBox = gobjects(maxEvt,1);
S.selectedLine = gobjects(maxEvt,1);
S.selectedText = gobjects(maxEvt,1);
S.boxLabel = gobjects(maxEvt,1);
S.finished = false;

S.fullXLim = [t(1) t(end)];
S.minWindowSec = max(5/Fs, 0.02);
S.zoomFactorIn = 0.4;
S.zoomFactorOut = 2.5;
S.rightClickHalfWidthSec = 20;

% Figure
fig = figure( ...
    'Name', 'Annotate switch events', ...
    'Color', 'w', ...
    'Units', 'normalized', ...
    'Position', [0.08 0.08 0.84 0.82], ...
    'NumberTitle', 'off', ...
    'WindowButtonDownFcn', @onMouseClickSafe, ...
    'WindowButtonMotionFcn', @onMouseMoveSafe, ...
    'WindowScrollWheelFcn', @onScrollSafe, ...
    'CloseRequestFcn', @onCloseSafe);

S.fig = fig;

ax = axes('Parent', fig, ...
    'Units', 'normalized', ...
    'Position', [0.07 0.12 0.60 0.80]);

S.ax = ax;

plot(ax, t, switchSig, 'k-', 'LineWidth', 1);
hold(ax, 'on');
ylim(ax, [0 1.5]);
xlim(ax, S.fullXLim);
xlabel(ax, 'Time (s)');
ylabel(ax, 'switchSig');
title(ax, 'Left click: assign event | Right click: x\pm20 s | Wheel: zoom x-axis');

% Moving cursor line
S.cursorLine = xline(ax, t(1), '--', 'LineWidth', 1, 'Color', [0.2 0.2 1]);
S.cursorLine.Visible = 'off';
S.cursorLine.HitTest = 'off';
S.cursorLine.PickableParts = 'none';

% Instruction text
uicontrol(fig, 'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.71 0.92 0.25 0.05], ...
    'String', 'always put start and end', ...
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'ForegroundColor', [0.8 0 0], ...
    'BackgroundColor', 'w');

uicontrol(fig, 'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.71 0.90 0.25 0.04], ...
    'String', 'Click a box to edit that event.', ...
    'FontSize', 10, ...
    'BackgroundColor', 'w');

uicontrol(fig, 'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.71 0.86 0.25 0.05], ...
    'String', 'Mouse wheel: zoom x-axis around cursor position', ...
    'FontSize', 10, ...
    'BackgroundColor', 'w');

uicontrol(fig, 'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.71 0.82 0.25 0.05], ...
    'String', 'Right click: show x from -20 to +20 s around cursor', ...
    'FontSize', 10, ...
    'BackgroundColor', 'w');

uicontrol(fig, 'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.71 0.78 0.25 0.05], ...
    'String', 'Recommend names: pair1start, pair1end, pair2start, pair2end, sleep1start, sleep1end', ...
    'FontSize', 10, ...
    'BackgroundColor', 'w');

% Event boxes
y0 = 0.75;
dy = 0.065;

for k = 1:maxEvt
    S.boxLabel(k) = uicontrol(fig, 'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.71 y0-(k-1)*dy 0.05 0.04], ...
        'String', sprintf('%d', k), ...
        'FontSize', 10, ...
        'BackgroundColor', 'w', ...
        'ForegroundColor', [0 0 0]);

    S.editBox(k) = uicontrol(fig, 'Style', 'edit', ...
        'Units', 'normalized', ...
        'Position', [0.77 y0-(k-1)*dy 0.18 0.045], ...
        'String', '', ...
        'FontSize', 10, ...
        'BackgroundColor', 'w', ...
        'UserData', k, ...
        'Callback', @onEditBoxClickedSafe, ...
        'ButtonDownFcn', @onEditBoxClickedSafe);
end

% Buttons
uicontrol(fig, 'Style', 'pushbutton', ...
    'Units', 'normalized', ...
    'Position', [0.72 0.08 0.10 0.05], ...
    'String', 'Finish', ...
    'FontSize', 10, ...
    'Callback', @onFinishSafe);

uicontrol(fig, 'Style', 'pushbutton', ...
    'Units', 'normalized', ...
    'Position', [0.84 0.08 0.10 0.05], ...
    'String', 'Clear last', ...
    'FontSize', 10, ...
    'Callback', @onClearLastSafe);

uicontrol(fig, 'Style', 'pushbutton', ...
    'Units', 'normalized', ...
    'Position', [0.72 0.015 0.10 0.05], ...
    'String', 'Reset axis', ...
    'FontSize', 10, ...
    'Callback', @onResetXSafe);

uicontrol(fig, 'Style', 'pushbutton', ...
    'Units', 'normalized', ...
    'Position', [0.84 0.015 0.10 0.05], ...
    'String', 'Clear all', ...
    'FontSize', 10, ...
    'Callback', @onClearAllSafe);

uicontrol(fig, 'Style', 'pushbutton', ...
    'Units', 'normalized', ...
    'Position', [0.96 0.015 0.03 0.05], ...
    'String', 'X', ...
    'FontSize', 10, ...
    'Callback', @onCloseSafe);

guidata(fig, S);
updateActiveBoxHighlight(fig);
uicontrol(S.editBox(1));

uiwait(fig);

if ~ishandle(fig)
    eventTable = table(strings(0,1), zeros(0,1), ...
        'VariableNames', {'eventName','signalIdx'});
    return
end

S = guidata(fig);

used = find(~isnan(S.selectedIdx));
eventNames = strings(numel(used),1);
signalIdx = S.selectedIdx(used);

for i = 1:numel(used)
    k = used(i);
    nm = strtrim(get(S.editBox(k), 'String'));
    if isempty(nm)
        nm = sprintf('event%d', k);
    end
    eventNames(i) = string(nm);
end

eventTable = table(eventNames, signalIdx, ...
    'VariableNames', {'eventName','signalIdx'});

if ~isempty(eventTable)
    [~, ord] = sort(eventTable.signalIdx);
    eventTable = eventTable(ord,:);
end

disp(eventTable);
return

    function onMouseMoveSafe(src, evt)
        try
            onMouseMove(src, evt);
        catch ME
            warning('%s', getReport(ME, 'basic'));
        end
    end

    function onMouseClickSafe(src, evt)
        try
            onMouseClick(src, evt);
        catch ME
            warning('%s', getReport(ME, 'basic'));
        end
    end

    function onScrollSafe(src, evt)
        try
            onScroll(src, evt);
        catch ME
            warning('%s', getReport(ME, 'basic'));
        end
    end

    function onResetXSafe(src, evt)
        try
            onResetX(src, evt);
        catch ME
            warning('%s', getReport(ME, 'basic'));
        end
    end

    function onEditBoxClickedSafe(src, evt)
        try
            onEditBoxClicked(src, evt);
        catch ME
            warning('%s', getReport(ME, 'basic'));
        end
    end

    function onClearLastSafe(src, evt)
        try
            onClearLast(src, evt);
        catch ME
            warning('%s', getReport(ME, 'basic'));
        end
    end

    function onClearAllSafe(src, evt)
        try
            onClearAll(src, evt);
        catch ME
            warning('%s', getReport(ME, 'basic'));
        end
    end

    function onFinishSafe(src, evt)
        try
            onFinish(src, evt);
        catch ME
            warning('%s', getReport(ME, 'basic'));
        end
    end

    function onCloseSafe(src, evt)
        try
            onClose(src, evt);
        catch ME
            warning('%s', getReport(ME, 'basic'));
            figLocal = ancestor(src, 'figure');
            if isempty(figLocal)
                figLocal = src;
            end
            if ishandle(figLocal)
                delete(figLocal);
            end
        end
    end

    function onMouseMove(src, ~)
        if ~ishandle(src)
            return
        end
        S = guidata(src);
        if isempty(S) || ~ishandle(S.ax)
            return
        end

        cp = get(S.ax, 'CurrentPoint');
        xNow = cp(1,1);
        yNow = cp(1,2);
        xl = xlim(S.ax);
        yl = ylim(S.ax);

        if xNow >= xl(1) && xNow <= xl(2) && yNow >= yl(1) && yNow <= yl(2)
            S.cursorLine.Value = xNow;
            S.cursorLine.Visible = 'on';
        else
            S.cursorLine.Visible = 'off';
        end
    end

    function onScroll(src, evt)
        if ~ishandle(src)
            return
        end
        S = guidata(src);
        if isempty(S) || ~ishandle(S.ax) || S.finished
            return
        end

        cp = get(S.ax, 'CurrentPoint');
        xNow = cp(1,1);
        yNow = cp(1,2);

        xl = xlim(S.ax);
        yl = ylim(S.ax);

        if xNow < xl(1) || xNow > xl(2) || yNow < yl(1) || yNow > yl(2)
            return
        end

        curWidth = diff(xl);

        if evt.VerticalScrollCount > 0
            newWidth = curWidth * S.zoomFactorOut;
        else
            newWidth = curWidth * S.zoomFactorIn;
        end

        fullWidth = diff(S.fullXLim);
        newWidth = min(newWidth, fullWidth);
        newWidth = max(newWidth, S.minWindowSec);

        frac = (xNow - xl(1)) / curWidth;
        if ~isfinite(frac)
            frac = 0.5;
        end

        newLeft = xNow - frac * newWidth;
        newRight = newLeft + newWidth;

        if newLeft < S.fullXLim(1)
            newLeft = S.fullXLim(1);
            newRight = newLeft + newWidth;
        end
        if newRight > S.fullXLim(2)
            newRight = S.fullXLim(2);
            newLeft = newRight - newWidth;
        end

        if newLeft < S.fullXLim(1)
            newLeft = S.fullXLim(1);
        end
        if newRight > S.fullXLim(2)
            newRight = S.fullXLim(2);
        end

        xlim(S.ax, [newLeft newRight]);
        S.cursorLine.Value = xNow;
        S.cursorLine.Visible = 'on';
    end

    function onMouseClick(src, ~)
        if ~ishandle(src)
            return
        end
        S = guidata(src);

        if S.finished
            return
        end

        cp = get(S.ax, 'CurrentPoint');
        xClick = cp(1,1);
        yClick = cp(1,2);

        xl = xlim(S.ax);
        yl = ylim(S.ax);
        if xClick < xl(1) || xClick > xl(2) || yClick < yl(1) || yClick > yl(2)
            return
        end

        clickType = get(src, 'SelectionType');

        if strcmp(clickType, 'alt')
            x1 = xClick - S.rightClickHalfWidthSec;
            x2 = xClick + S.rightClickHalfWidthSec;

            x1 = max(x1, S.fullXLim(1));
            x2 = min(x2, S.fullXLim(2));

            if (x2 - x1) < 2*S.rightClickHalfWidthSec
                if x1 <= S.fullXLim(1)
                    x2 = min(S.fullXLim(2), x1 + 2*S.rightClickHalfWidthSec);
                elseif x2 >= S.fullXLim(2)
                    x1 = max(S.fullXLim(1), x2 - 2*S.rightClickHalfWidthSec);
                end
            end

            xlim(S.ax, [x1 x2]);
            S.cursorLine.Value = xClick;
            S.cursorLine.Visible = 'on';
            return
        end

        k = S.activeBox;
        if k < 1 || k > S.maxEvt
            return
        end

        thisIdx = round(xClick * S.Fs) + 1;
        thisIdx = max(1, min(S.nSamp, thisIdx));
        thisTime = (thisIdx - 1) / S.Fs;

        if isgraphics(S.selectedLine(k))
            delete(S.selectedLine(k));
        end
        if isgraphics(S.selectedText(k))
            delete(S.selectedText(k));
        end

        S.selectedIdx(k) = thisIdx;
        S.selectedLine(k) = xline(S.ax, thisTime, 'b-', 'LineWidth', 2);
        S.selectedLine(k).HitTest = 'off';
        S.selectedLine(k).PickableParts = 'none';

        S.selectedText(k) = text(S.ax, thisTime, 1.05, sprintf('event%d', k), ...
            'Color', 'b', ...
            'Rotation', 90, ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'left', ...
            'FontSize', 10, ...
            'FontWeight', 'bold', ...
            'HitTest', 'off', ...
            'PickableParts', 'none');

        if k < S.maxEvt
            S.activeBox = k + 1;
        end

        guidata(src, S);
        updateActiveBoxHighlight(src);

        S = guidata(src);
        if S.activeBox <= S.maxEvt
            uicontrol(S.editBox(S.activeBox));
        end
    end

    function onResetX(src, ~)
        figLocal = ancestor(src, 'figure');
        if isempty(figLocal)
            figLocal = src;
        end
        if ~ishandle(figLocal)
            return
        end
        S = guidata(figLocal);
        if isempty(S) || ~ishandle(S.ax)
            return
        end
        xlim(S.ax, S.fullXLim);
    end

    function onEditBoxClicked(src, ~)
        figLocal = ancestor(src, 'figure');
        S = guidata(figLocal);
        if S.finished
            return
        end
        k = get(src, 'UserData');
        S.activeBox = k;
        guidata(figLocal, S);
        updateActiveBoxHighlight(figLocal);
        uicontrol(S.editBox(k));
    end

    function onClearLast(src, ~)
        figLocal = ancestor(src, 'figure');
        S = guidata(figLocal);
        if S.finished
            return
        end

        used = find(~isnan(S.selectedIdx));
        if isempty(used)
            return
        end

        k = used(end);

        if isgraphics(S.selectedLine(k))
            delete(S.selectedLine(k));
        end
        if isgraphics(S.selectedText(k))
            delete(S.selectedText(k));
        end

        S.selectedIdx(k) = NaN;
        set(S.editBox(k), 'String', '');
        S.activeBox = k;

        guidata(figLocal, S);
        updateActiveBoxHighlight(figLocal);
        uicontrol(S.editBox(k));
    end

    function onClearAll(src, ~)
        figLocal = ancestor(src, 'figure');
        S = guidata(figLocal);
        if S.finished
            return
        end

        for k = 1:S.maxEvt
            if isgraphics(S.selectedLine(k))
                delete(S.selectedLine(k));
            end
            if isgraphics(S.selectedText(k))
                delete(S.selectedText(k));
            end
            S.selectedIdx(k) = NaN;
            set(S.editBox(k), 'String', '');
        end

        S.activeBox = 1;
        guidata(figLocal, S);
        updateActiveBoxHighlight(figLocal);
        uicontrol(S.editBox(1));
    end

    function onFinish(src, ~)
        figLocal = ancestor(src, 'figure');
        if isempty(figLocal)
            figLocal = src;
        end
        if ~ishandle(figLocal)
            return
        end
        S = guidata(figLocal);
        S.finished = true;
        guidata(figLocal, S);

        set(figLocal, ...
            'WindowButtonDownFcn', '', ...
            'WindowButtonMotionFcn', '', ...
            'WindowScrollWheelFcn', '');

        uiresume(figLocal);
    end

    function onClose(src, ~)
        figLocal = ancestor(src, 'figure');
        if isempty(figLocal)
            figLocal = src;
        end
        if ~ishandle(figLocal)
            return
        end
        delete(figLocal);
    end

    function updateActiveBoxHighlight(figLocal)
        if ~ishandle(figLocal)
            return
        end
        S = guidata(figLocal);
        for kk = 1:S.maxEvt
            if kk == S.activeBox && ~S.finished
                set(S.editBox(kk), 'BackgroundColor', [1 1 0.85]);
                set(S.boxLabel(kk), 'ForegroundColor', [0 0 1], 'FontWeight', 'bold');
            else
                set(S.editBox(kk), 'BackgroundColor', 'w');
                set(S.boxLabel(kk), 'ForegroundColor', [0 0 0], 'FontWeight', 'normal');
            end
        end
    end
end