import tensorflow as tf
import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score

class PhagePredictor:
    """Phage Classifier - For Loading and Using Trained Models"""
    
    def __init__(self, model_path, seq_length=20, nb_classes=2):
        """
        Initialize the predictor

        Parameters:
        model_path: Path to the trained model file (.h5)
        seq_length: Sequence length (default is 20, corresponding to 2 classes * 10 samples)
        nb_classes: Number of classes (default is 2)
        """

        self.model_path = model_path
        self.seq_length = seq_length
        self.nb_classes = nb_classes
        
        # Load the model
        self.model = self.load_model()
        
        # Print model information
        self.print_model_info()
    
    def load_model(self):
        """Load the trained model"""
        if not os.path.exists(self.model_path):
            raise FileNotFoundError(f"Model file does not exist: {self.model_path}")
        
        print(f"Loading model: {self.model_path}")
        model = tf.keras.models.load_model(self.model_path)
        print("Model loaded successfully!")
        return model
    
    def print_model_info(self):
        """Print model information"""
        print("\n=== Model Information ===")
        print(f"Input shape: {self.model.input_shape}")
        print(f"Output shape: {self.model.output_shape}")
        print(f"Number of layers: {len(self.model.layers)}")
        
        # Display model architecture summary
        print("\nModel architecture:")
        self.model.summary()
    
    def prepare_input_data(self, npy_file_path):
        """
        Prepare input data

        Parameters:
        npy_file_path: Path to the npy feature file

        Returns:
        Processed data with shape (n_samples, seq_length, feature_dim)
        """
        # Load npy file
        print(f"Loading feature file: {npy_file_path}")
        features = np.load(npy_file_path)
        
        # Check feature dimensions
        print(f"Original data shape: {features.shape}")
        
        # Determine feature dimensions
        if len(features.shape) == 1:
            # If it's a 1D array, reshape it to 2D (1, n_features)
            features = features.reshape(1, -1)
        elif len(features.shape) > 2:
            raise ValueError(f"Input data dimension not supported: {features.shape}")
        
        n_samples, feature_dim = features.shape
        print(f"Number of samples: {n_samples}, feature dimension: {feature_dim}")
        
        # To be compatible with the sequence model used during training, create sequence data
        # Repeat each sample seq_length times to form a sequence
        prepared_data = []
        
        for i in range(n_samples):
            # Create a sequence by repeating the same sample seq_length times
            sequence = np.repeat(features[i:i+1], self.seq_length, axis=0)
            prepared_data.append(sequence)
        
        prepared_data = np.array(prepared_data)
        print(f"Processed data shape: {prepared_data.shape}")
        
        return prepared_data, features
    
    def predict(self, npy_file_path, threshold=0.5):
        """
        Predict using the input npy file

        Parameters:
        npy_file_path: Path to the npy feature file
        threshold: Classification threshold (default 0.5)

        Returns:
        Prediction result dictionary containing probabilities, classes, and detailed information
        """
        # Prepare data
        prepared_data, original_features = self.prepare_input_data(npy_file_path)
        
        # Make predictions
        print("\nMaking predictions..")
        predictions = self.model.predict(prepared_data, verbose=1)
        
        # Process prediction results
        # Model output shape: (n_samples, seq_length, nb_classes)
        # We take the prediction at the last time step for each sequence
        last_step_predictions = predictions[:, -1, :]
        
        # Compute probabilities
        probabilities = tf.nn.softmax(last_step_predictions).numpy()
        
        # Extract probability for lytic phage (class 1)
        lytic_probabilities = probabilities[:, 1]
        
        # Determine classes based on threshold
        classes = (lytic_probabilities >= threshold).astype(int)
        
        # Create detailed prediction results
        results = []
        for i in range(len(lytic_probabilities)):
            lytic_prob = lytic_probabilities[i]
            class_label = classes[i]
            
            # Determine class name
            if class_label == 0:
                class_name = "Temperate phage"
                confidence = 1 - lytic_prob  # Confidence for temperate phage
            else:
                class_name = "Lytic phage"
                confidence = lytic_prob  # Confidence for lytic phage
            
            result = {
                'Sample index': i,
                'Lytic probability': float(lytic_prob),
                'Temperate probability': float(1 - lytic_prob),
                'Predicted class': class_label,
                'Class name': class_name,
                'Confidence': float(confidence),
                'Threshold': threshold
            }
            results.append(result)
        
        # Summary statistics
        total_samples = len(results)
        temperate_count = sum(1 for r in results if r['Predicted class'] == 0)
        lytic_count = sum(1 for r in results if r['Predicted class'] == 1)
        
        summary = {
            'Total samples': total_samples,
            'Number of temperate phages': temperate_count,
            'Number of lytic phages': lytic_count,
            'Proportion of temperate phages': temperate_count / total_samples if total_samples > 0 else 0,
            'Proportion of lytic phages': lytic_count / total_samples if total_samples > 0 else 0,
            'Threshold used': threshold,
            'Detailed results': results
        }
        
        return summary, probabilities, original_features
    
    def save_predictions(self, summary, output_dir='./predictions'):
        """Save prediction results to file"""
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = np.datetime64('now').astype(str).replace(':', '-').replace(' ', '_')
        output_file = os.path.join(output_dir, f'predictions_{timestamp}.txt')
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("=== Phage Classification Prediction Results ===\n\n")
            f.write(f"Model path: {self.model_path}\n")
            f.write(f"Prediction timestamp: {timestamp}\n")
            f.write(f"Total samples: {summary['Total samples']}\n")
            f.write(f"Number of temperate phages: {summary['Number of temperate phages']} ({summary['Proportion of temperate phages']*100:.2f}%)\n")
            f.write(f"Number of lytic phages: {summary['Number of lytic phages']} ({summary['Proportion of lytic phages']*100:.2f}%)\n")
            f.write(f"Classification threshold: {summary['Threshold used']}\n")
            
            f.write("\n=== Detailed Prediction Results ===\n")
            for result in summary['Detailed results']:
                f.write(f"\nSample {result['Sample index']}:\n")
                f.write(f"  Predicted class: {result['Class name']}\n")
                f.write(f"  Lytic probability: {result['Lytic probability']:.4f}\n")
                f.write(f"  Temperate probability: {result['Temperate probability']:.4f}\n")
                f.write(f"  Confidence: {result['Confidence']:.4f}\n")
        
        print(f"Prediction results have been saved to: {output_file}")
        
        return output_file
    
    def plot_probability_distribution(self, probabilities, output_dir='./predictions'):
        """Plot probability distribution"""
        # Extract lytic probabilities
        lytic_probs = probabilities[:, 1]
        
        plt.figure(figsize=(10, 6))
        
        # Create histogram
        plt.hist(lytic_probs, bins=50, alpha=0.7, color='steelblue', edgecolor='black')
        
        # Add threshold line
        plt.axvline(x=0.5, color='red', linestyle='--', linewidth=2, label='Classification threshold (0.5)')
        
        # Add statistics line
        mean_prob = np.mean(lytic_probs)
        std_prob = np.std(lytic_probs)
        plt.axvline(x=mean_prob, color='green', linestyle='-', linewidth=2, label=f'Mean ({mean_prob:.3f})')
        
        plt.xlabel('Predicted probability of lytic phage', fontsize=12)
        plt.ylabel('Number of samples', fontsize=12)
        plt.title('Predicted Probability Distribution', fontsize=14, fontweight='bold')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Add statistics text box
        stats_text = f'Number of samples: {len(lytic_probs)}\n' \
                    f'Mean: {mean_prob:.3f}\n' \
                    f'Standard deviation: {std_prob:.3f}\n' \
                    f'Min: {np.min(lytic_probs):.3f}\n' \
                    f'Max: {np.max(lytic_probs):.3f}'
        
        plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes,
                verticalalignment='top', horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # Save plot
        os.makedirs(output_dir, exist_ok=True)
        timestamp = np.datetime64('now').astype(str).replace(':', '-').replace(' ', '_')
        plot_file = os.path.join(output_dir, f'probability_distribution_{timestamp}.png')
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"Probability distribution plot has been saved to: {plot_file}")
        
        return plot_file
    
    def batch_predict(self, npy_files, output_dir='./batch_predictions'):
        """
        Batch predict multiple npy files

        Parameters:
        npy_files: List of npy file paths
        output_dir: Output directory

        Returns:
        Summary of prediction results for all files
        """
        os.makedirs(output_dir, exist_ok=True)
        
        all_results = []
        all_probabilities = []
        
        for i, npy_file in enumerate(npy_files):
            print(f"\nProcessing file {i+1}/{len(npy_files)}: {npy_file}")
            
            try:
                summary, probabilities, _ = self.predict(npy_file)
                
                # Save prediction results for a single file
                file_results = {
                    'File name': os.path.basename(npy_file),
                    'File path': npy_file,
                    'Number of samples': summary['Total samples'],
                    'Number of temperate phages': summary['Number of temperate phages'],
                    'Number of lytic phages': summary['Number of lytic phages'],
                    'Proportion of temperate phages': summary['Proportion of temperate phages'],
                    'Proportion of lytic phages': summary['Proportion of lytic phages']
                }
                all_results.append(file_results)
                all_probabilities.extend(probabilities[:, 1].tolist())
                
                # Save detailed results for the single file
                single_file_dir = os.path.join(output_dir, 'single_files')
                os.makedirs(single_file_dir, exist_ok=True)
                self.save_predictions(
                    summary, 
                    os.path.join(single_file_dir, f"file_{i+1}")
                )
                
            except Exception as e:
                print(f"Error processing file {npy_file}: {str(e)}")
                continue
        
        # Save batch prediction summary
        if all_results:
            self.save_batch_summary(all_results, all_probabilities, output_dir)
        
        return all_results
    
    def save_batch_summary(self, results, probabilities, output_dir):
        """Save batch prediction summary"""
        summary_file = os.path.join(output_dir, 'batch_summary.txt')
        
        total_samples = sum(r['Number of samples'] for r in results)
        total_temperate = sum(r['Number of temperate phages'] for r in results)
        total_lytic = sum(r['Number of lytic phages'] for r in results)
        
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("=== Batch Prediction Summary ===\n\n")
            f.write(f"Model path: {self.model_path}\n")
            f.write(f"Number of files: {len(results)}\n")
            f.write(f"Total samples: {total_samples}\n")
            f.write(f"Total temperate phages: {total_temperate} ({total_temperate/total_samples*100:.2f}%)\n")
            f.write(f"Total lytic phages: {total_lytic} ({total_lytic/total_samples*100:.2f}%)\n")
            
            f.write("\n=== Detailed information for each file ===\n")
            for i, result in enumerate(results):
                f.write(f"\nFile {i+1}: {result['File name']}\n")
                f.write(f"  Number of samples: {result['Number of samples']}\n")
                f.write(f"  Number of temperate phages: {result['Number of temperate phages']} ({result['Proportion of temperate phages']*100:.2f}%)\n")
                f.write(f"  Number of lytic phages: {result['Number of lytic phages']} ({result['Proportion of lytic phages']*100:.2f}%)\n")
        
        print(f"Batch prediction summary has been saved to: {summary_file}")
        
        # Plot probability distribution for all samples
        if probabilities:
            plt.figure(figsize=(10, 6))
            plt.hist(probabilities, bins=50, alpha=0.7, color='steelblue', edgecolor='black')
            plt.axvline(x=0.5, color='red', linestyle='--', linewidth=2, label='Classification threshold (0.5)')
            plt.xlabel('Predicted probability of lytic phage')
            plt.ylabel('Number of samples')
            plt.title('Batch Prediction Probability Distribution')
            plt.legend()
            plt.grid(True, alpha=0.3)
            
            plot_file = os.path.join(output_dir, 'batch_probability_distribution.png')
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.show()
            
            print(f"Batch prediction probability distribution plot has been saved to: {plot_file}")


# Usage example
if __name__ == "__main__":
    # ========== Configuration Parameters ==========
    MODEL_PATH = "./final_model.h5"  # Change to your model path
    NPY_FILE_PATH = "./your_new_data.npy"  # Change to your npy file path
    
    # ========== Create Predictor ==========
    predictor = PhagePredictor(
        model_path=MODEL_PATH,
        seq_length=20,  # Adjust according to your training settings
        nb_classes=2    # Binary classification
    )
    
    # ========== Make Predictions ==========
    # Single file prediction
    print("\n" + "="*50)
    print("Starting single file prediction")
    print("="*50)
    
    summary, probabilities, features = predictor.predict(NPY_FILE_PATH, threshold=0.5)
    
    # Print results
    print("\n" + "="*50)
    print("Prediction Summary")
    print("="*50)
    print(f"Total samples: {summary['Total samples']}")
    print(f"Number of temperate phages: {summary['Number of temperate phages']} ({summary['Proportion of temperate phages']*100:.2f}%)")
    print(f"Number of lytic phages: {summary['Number of lytic phages']} ({summary['Proportion of lytic phages']*100:.2f}%)")
    print(f"Classification threshold: {summary['Threshold used']}")
    
    # Save prediction results
    output_file = predictor.save_predictions(summary)
    
    # Plot probability distribution
    plot_file = predictor.plot_probability_distribution(probabilities)
    
    # ========== Batch Prediction Example ==========
    # If you have multiple npy files, you can use the following code for batch prediction
    """
    npy_files = [
        "./data/file1.npy",
        "./data/file2.npy",
        "./data/file3.npy"
    ]
    
    batch_results = predictor.batch_predict(npy_files, output_dir="./batch_results")
    """
    
    print("\n" + "="*50)
    print("Prediction completed!")
    print("="*50)