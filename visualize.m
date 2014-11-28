function visualize(filename)
	file = fopen(filename);
	gridSize = fscanf(file, '%d %d');
	gridSize = gridSize(1:2);
	entries = fscanf(file, '%f %f %f %f %f', [5,Inf]);
	x = reshape(entries(1,:),[gridSize(2),gridSize(1)]).';
	y = reshape(entries(2,:),[gridSize(2),gridSize(1)]).';
	u = reshape(entries(3,:),[gridSize(2),gridSize(1)]).';
	v = reshape(entries(4,:),[gridSize(2),gridSize(1)]).';
	p = reshape(entries(5,:),[gridSize(2),gridSize(1)]).';

	subplot(2,2,1);
	surf(x,y,u,'EdgeColor','None');
	title('u');
	xlabel('x');
	xlabel('y');

	subplot(2,2,2);
	surf(x,y,v,'EdgeColor','None');
	title('v');
	xlabel('x');
	xlabel('y');

	subplot(2,2,3);
	surf(x,y,p,'EdgeColor','None');
	title('p');
	xlabel('x');
	xlabel('y');
    
    subplot(2,2,4);
    startx = rand([1,floor(gridSize(1)*1.7)]);
    streamline(x,y,u,v,startx,rand([1,size(startx)]));
    title('streamlines');
    
   % figure(2);
   % step = 8;
	%quiver(x(1:step:end),y(1:step:end),u(1:step:end),v(1:step:end));
	%title('p');
	%xlabel('x');
	%xlabel('y');
end
