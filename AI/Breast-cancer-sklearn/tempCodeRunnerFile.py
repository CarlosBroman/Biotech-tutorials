sample = np.array([[8,1,1,8,2,1,9,7,1]]) # Characteristics based on the data, 
                                            # change these numbers for different predictions
result = loaded_model.predict(sample)
print(classes[int(result)])