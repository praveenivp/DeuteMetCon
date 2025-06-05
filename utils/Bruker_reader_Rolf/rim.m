function rim(data)
% Calls ImageWindow
image = ImageWindow(data);
if numel(image.widgets) == 0
    delete(image);
    return;
end
name = inputname(1);
if numel(name)>0
    image.widgets.base.Name=name;
else
    s = size(data);
    name = [num2str(s(1)),' x ',num2str(s(2))];
    if numel(s) > 2
        name = [name, ' x ',num2str(s(3))];
    end
    image.widgets.base.Name=name;
end

