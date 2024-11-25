# Example training data
data(iris)
training_data <- iris[, -5]
training_labels <- iris$Species

# Train the SVM
svm_model <- svm(training_data, training_labels, kernel = "radial", probability = TRUE)

# New test data (ensure it matches training features)
test_data <- iris[, -5]  # Matching column structure

# Predict with proper feature alignment
y_pred <- predict(svm_model, newdata = test_data, decision.values = TRUE, probability = TRUE)