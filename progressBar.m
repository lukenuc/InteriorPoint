function [] = progressBar(percentComplete)

    % Persistent variables to store the current fun fact and the last displayed percentage
    persistent currentFunFact;
    persistent lastPercent;

    if isempty(currentFunFact)
        currentFunFact = '';
    end

    if isempty(lastPercent)
        lastPercent = -20;
    end

    % Fun facts
    funFacts = { ...
        'A "jiffy" is an actual unit of time. It''s 1/100th of a second.', ...
        'There are more stars in the universe than grains of sand on all the Earth''s beaches.', ...
        'The world consumes about 2.25 billion cups of coffee every day. That''s a lot of energy!', ...
        'A prime number is a natural number greater than 1 that cannot be formed by multiplying two smaller natural numbers.', ...
        'Your brain uses about 20% of the total oxygen and blood in your body, even though it makes up only about 2% of your body weight.', ...
        'Coffee beans aren''t actually beans. They are the seeds of coffee cherries.', ...
        'The number "pi" is not just a tasty dessert; it''s the ratio of a circle''s circumference to its diameter.', ...
        'Honey never spoils. Archaeologists have found pots of honey in ancient Egyptian tombs that are over 3,000 years old and still edible.', ...
        'Beethoven counted exactly 60 coffee beans per cup when he brewed his coffee.'};

    % Calculate the number of '=' signs to display
    numEquals = floor(percentComplete / 2);

    % Check if it's time to display a new fun fact
    if percentComplete >= lastPercent + 20
        currentFunFact = funFacts{randi(length(funFacts))};
        lastPercent = floor(percentComplete / 20) * 20;
    end

    % Print the progress bar with the current fun fact
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    fprintf('[%s%s] %3.0f%% \n%s\n', repmat('=', 1, numEquals), repmat(' ', 1, 50-numEquals), percentComplete, currentFunFact);

end
