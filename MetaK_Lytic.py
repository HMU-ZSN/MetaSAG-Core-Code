import tensorflow as tf
from tensorflow.keras import layers, models
import argparse
import os
import numpy as np
import datetime
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
import seaborn as sns

# Set the Chinese font and graphic style
plt.rcParams['font.sans-serif'] = ['SimHei']  # It is used to display Chinese tags normally
plt.rcParams['axes.unicode_minus'] = False  # It is used to display the minus sign normally
sns.set_style("whitegrid")


class PhageLSTMModel:

    def __init__(self, input_size, hidden_size, num_layers, num_classes,
                 seq_length, model_type="LSTM"):
        self.input_size = input_size
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        self.num_classes = num_classes
        self.seq_length = seq_length
        self.model_type = model_type

        self.model = self.build_model()

    def build_model(self):
        # Build the model using the Keras Functional API
        inputs = tf.keras.Input(shape=(self.seq_length, self.input_size), name="image")

        # Build the RNN layer
        if self.model_type == "LSTM":
            rnn_layer = layers.LSTM(self.hidden_size, return_sequences=True)
        else:  # GRU
            rnn_layer = layers.GRU(self.hidden_size, return_sequences=True)

        x = inputs
        # Multi-layer RNN
        for _ in range(self.num_layers - 1):
            if self.model_type == "LSTM":
                x = layers.LSTM(self.hidden_size, return_sequences=True)(x)
            else:  # GRU
                x = layers.GRU(self.hidden_size, return_sequences=True)(x)

        # The last layer is RNN
        x = rnn_layer(x)

        # output layer
        outputs = layers.Dense(self.num_classes, activation=None)(x)

        model = tf.keras.Model(inputs=inputs, outputs=outputs)
        return model


class PhageFeatureGenerator:
    """Phage feature data generator"""

    def __init__(self, temperate_file, lytic_file, nb_classes=2, nb_samples_per_class=10, train_ratio=0.8):
        self.nb_classes = nb_classes
        self.nb_samples_per_class = nb_samples_per_class

        # Load the feature file
        print(f"Loading temperate phage features from: {temperate_file}")
        temperate_features = np.load(temperate_file)
        print(f"Loading lytic phage features from: {lytic_file}")
        lytic_features = np.load(lytic_file)

        print(f"Temperate phage features shape: {temperate_features.shape}")
        print(f"Lytic phage features shape: {lytic_features.shape}")

        # Ensure consistency in feature dimensions
        assert temperate_features.shape[1] == lytic_features.shape[1], "Feature dimensions must match"
        self.feature_dim = temperate_features.shape[1]
        print(f"Feature dimension: {self.feature_dim}")

        # Creating tags
        temperate_labels = np.zeros(len(temperate_features))  # The label of the mild bacteriophage is 0
        lytic_labels = np.ones(len(lytic_features))  # The label for virulent bacteriophages is 1

        # Merge features and labels
        all_features = np.vstack([temperate_features, lytic_features])
        all_labels = np.hstack([temperate_labels, lytic_labels])

        # Disrupt the data
        indices = np.random.permutation(len(all_features))
        all_features = all_features[indices]
        all_labels = all_labels[indices]

        # Split the training set and the test set
        split_idx = int(train_ratio * len(all_features))
        self.train_features = all_features[:split_idx]
        self.train_labels = all_labels[:split_idx]
        self.test_features = all_features[split_idx:]
        self.test_labels = all_labels[split_idx:]

        print(f"Training samples: {len(self.train_features)}")
        print(f"Test samples: {len(self.test_features)}")
        print(f"Temperate samples in training: {np.sum(self.train_labels == 0)}")
        print(f"Lytic samples in training: {np.sum(self.train_labels == 1)}")

    def sample(self, mode="train", batch_size=16):
        """Generate batch data"""
        if mode == "train":
            features = self.train_features
            labels = self.train_labels
        else:
            features = self.test_features
            labels = self.test_labels

        # Randomly select batch samples
        indices = np.random.choice(len(features), batch_size, replace=False)
        batch_features = features[indices]
        batch_labels = labels[indices]

        # To be compatible with RNN, we need to create sequence data
        # Each sequence contains samples of nb_classes * nb_samples_per_class
        seq_length = self.nb_classes * self.nb_samples_per_class
        sequence_features = []
        sequence_labels = []

        for i in range(batch_size):
            # Create a sequence for each batch of samples
            seq_feat = []
            seq_lab = []

            # Randomly select samples to construct the sequence
            for j in range(seq_length):
                idx = np.random.randint(len(features))
                seq_feat.append(features[idx])
                seq_lab.append(labels[idx])

            sequence_features.append(seq_feat)
            sequence_labels.append(seq_lab)

        return np.array(sequence_features), np.array(sequence_labels)

    def get_all_test_data(self):
        """Obtain all test data for the final evaluation"""
        # To be compatible with RNN, we need to create sequence data
        seq_length = self.nb_classes * self.nb_samples_per_class
        test_features = []
        test_labels = []

        # Create sufficient sequences to cover all test samples
        num_sequences = len(self.test_features) // seq_length + 1

        for i in range(num_sequences):
            seq_feat = []
            seq_lab = []

            for j in range(seq_length):
                idx = (i * seq_length + j) % len(self.test_features)
                seq_feat.append(self.test_features[idx])
                seq_lab.append(self.test_labels[idx])

            test_features.append(seq_feat)
            test_labels.append(seq_lab)

        return np.array(test_features), np.array(test_labels)


def build_argparser():
    parser = argparse.ArgumentParser()

    parser.add_argument('--mode', default="train")
    parser.add_argument('--batch-size',
                        dest='_batch_size', help='Batch size (default: %(default)s)',
                        type=int, default=16)
    parser.add_argument('--num-classes',
                        dest='_nb_classes', help='Number of classes (default: %(default)s)',
                        type=int, default=2)
    parser.add_argument('--num-samples',
                        dest='_nb_samples_per_class', help='Number of samples in each episode (default: %(default)s)',
                        type=int, default=10)
    parser.add_argument('--hidden-size',
                        dest='_hidden_size', help='Number of hidden units in RNN (default: %(default)s)',
                        type=int, default=200)
    parser.add_argument('--num_layers',
                        dest='_num_layers', help='Number of RNN layers (default: %(default)s)',
                        type=int, default=1)
    parser.add_argument('--learning-rate',
                        dest='_learning_rate', help='Learning Rate (default: %(default)s)',
                        type=float, default=1e-3)
    parser.add_argument('--iterations',
                        dest='_iterations', help='Number of iterations for training (default: %(default)s)',
                        type=int, default=100000)
    parser.add_argument('--temperate-file', default='./temperate_features.npy',
                        help='File containing temperate phage features')
    parser.add_argument('--lytic-file', default='./lytic_features.npy',
                        help='File containing lytic phage features')
    parser.add_argument('--train-ratio', default=0.8, type=float,
                        help='Ratio of data to use for training (default: %(default)s)')
    parser.add_argument('--save-dir', default='./ckpt/')
    parser.add_argument("--log-dir", default="./log/")
    parser.add_argument('--model', default="LSTM", help='LSTM or GRU')

    return parser


def metric_accuracy(args, labels, outputs):
    """Calculate the accuracy index"""
    seq_length = args._nb_classes * args._nb_samples_per_class
    outputs = np.argmax(outputs, axis=-1)
    correct = [0] * seq_length
    total = [0] * seq_length
    for i in range(np.shape(labels)[0]):
        label = labels[i]
        output = outputs[i]
        class_count = {}
        for j in range(seq_length):
            class_count[label[j]] = class_count.get(label[j], 0) + 1
            total[class_count[label[j]]] += 1
            if label[j] == output[j]:
                correct[class_count[label[j]]] += 1
    return [float(correct[i]) / total[i] if total[i] > 0. else 0. for i in range(1, args._nb_samples_per_class + 1)]


def calculate_auc_aupr(model, features, labels):
    """Calculate the UC-ROC and AUPR values"""
    
    predictions = model.predict(features, verbose=0)

    
    seq_length = predictions.shape[1]
    flat_predictions = predictions.reshape(-1, predictions.shape[-1])
    flat_labels = labels.reshape(-1)

    
    positive_probs = tf.nn.softmax(flat_predictions).numpy()[:, 1]

    
    fpr, tpr, _ = roc_curve(flat_labels, positive_probs)
    roc_auc = auc(fpr, tpr)

    
    precision, recall, _ = precision_recall_curve(flat_labels, positive_probs)
    aupr = average_precision_score(flat_labels, positive_probs)

    return fpr, tpr, roc_auc, precision, recall, aupr, positive_probs, flat_labels


def plot_metrics(fpr, tpr, roc_auc, precision, recall, aupr, save_path):
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    
    ax1.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {roc_auc:.4f})')
    ax1.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', label='Random classifier')
    ax1.set_xlim([0.0, 1.0])
    ax1.set_ylim([0.0, 1.05])
    ax1.set_xlabel('false positive rate (False Positive Rate)')
    ax1.set_ylabel('True Positive Rate (True Positive Rate)')
    ax1.set_title(' ROC')
    ax1.legend(loc="lower right")
    ax1.grid(True)

   
    ax2.plot(recall, precision, color='blue', lw=2, label=f'PR curve (AUPR = {aupr:.4f})')
    ax2.set_xlim([0.0, 1.0])
    ax2.set_ylim([0.0, 1.05])
    ax2.set_xlabel('Recall')
    ax2.set_ylabel('Precision')
    ax2.set_title('Precision-Recall')
    ax2.legend(loc="upper right")
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"The chart has been saved to: {save_path}")
    print(f"AUC-ROC: {roc_auc:.4f}")
    print(f"AUPR: {aupr:.4f}")


def train(model, data_generator, args):
    
    model.compile(
        optimizer=tf.keras.optimizers.Adam(learning_rate=args._learning_rate),
        loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True),
        metrics=['accuracy']
    )

   
    os.makedirs(os.path.join(args.save_dir, args.model), exist_ok=True)
    os.makedirs(args.log_dir, exist_ok=True)

   
    callbacks = [
        tf.keras.callbacks.ModelCheckpoint(
            os.path.join(args.save_dir, args.model, "model.h5"),
            save_best_only=False,  
            save_weights_only=False
        ),
        tf.keras.callbacks.CSVLogger(
            os.path.join(args.log_dir, f"{args.model}-{datetime.datetime.now().strftime('%m-%d-%H-%M')}.csv")
        )
    ]

    print(args)
    print("Training model...")

 
    best_auc = 0
    best_aupr = 0
    best_iteration = 0

    # Customize the training loop to match the original logic
    for ep in range(args._iterations):
        if ep % 100 == 0:
            # Evaluate on the test set
            test_features, test_labels = data_generator.sample("test", args._batch_size)
            test_loss, test_accuracy = model.evaluate(test_features, test_labels, verbose=0)

            
            test_outputs = model.predict(test_features, verbose=0)
            accuracy = metric_accuracy(args, test_labels, test_outputs)

            
            fpr, tpr, roc_auc, precision, recall, aupr, _, _ = calculate_auc_aupr(
                model, test_features, test_labels
            )

            # Update to the best performance
            if roc_auc > best_auc:
                best_auc = roc_auc
                best_aupr = aupr
                best_iteration = ep
                # Save the best model
                model.save(os.path.join(args.save_dir, args.model, "best_model.h5"))

            for accu in accuracy:
                print('%.4f' % accu, end='\t')
            print('%d\t%.4f\t%.4f\t%.4f' % (ep, test_loss, roc_auc, aupr))

        # training steps

        train_features, train_labels = data_generator.sample("train", args._batch_size)
        model.train_on_batch(train_features, train_labels)

    # Save the final model after the training is completed
    model.save(os.path.join(args.save_dir, args.model, "final_model.h5"))

    # Output the best performance
    print(f"\nTraining completed!")
    print(f"Optimal AUC-ROC: {best_auc:.4f} (in iteration {best_iteration})")
    print(f"Best AUPR: {best_aupr:.4f} (in iteration {best_iteration})")


def test(model, data_generator, args):
    
    print("Test the model performance...")

    
    test_features, test_labels = data_generator.get_all_test_data()

  
    test_loss, test_accuracy = model.evaluate(test_features, test_labels, verbose=0)

    
    fpr, tpr, roc_auc, precision, recall, aupr, probs, flat_labels = calculate_auc_aupr(
        model, test_features, test_labels
    )

    
    test_outputs = model.predict(test_features, verbose=0)
    accuracy = metric_accuracy(args, test_labels, test_outputs)

    
    print("\ntest result:")
    print("1st\t2nd\t3rd\t4th\t5th\t6th\t7th\t8th\t9th\t10th\tloss\tAUC-ROC\tAUPR")
    for accu in accuracy:
        print('%.4f' % accu, end='\t')
    print('%.4f\t%.4f\t%.4f' % (test_loss, roc_auc, aupr))

    
    plot_save_path = os.path.join(args.log_dir,
                                  f"{args.model}_metrics_{datetime.datetime.now().strftime('%m-%d-%H-%M')}.pdf")
    plot_metrics(fpr, tpr, roc_auc, precision, recall, aupr, plot_save_path)

    
    plot_probability_distribution(probs, flat_labels, args.log_dir, args.model)

    return roc_auc, aupr


def plot_probability_distribution(probs, labels, log_dir, model_name):
    
    plt.figure(figsize=(10, 6))

    
    temperate_probs = probs[labels == 0]  
    lytic_probs = probs[labels == 1]  

    
    bins = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
    bin_labels = ['0-0.2', '0.2-0.4', '0.4-0.6', '0.6-0.8', '0.8-1.0']

   
    print("\n=== Probability distribution statistics ===")

    # Statistics of mild phages
    temperate_counts, _ = np.histogram(temperate_probs, bins=bins)
    print("\nProbability distribution of lysogenic phages:")
    for i, (bin_label, count) in enumerate(zip(bin_labels, temperate_counts)):
        print(f"  {bin_label}: {count} samples ({count / len(temperate_probs) * 100:.2f}%)")

    # Statistics of virulent bacteriophages
    lytic_counts, _ = np.histogram(lytic_probs, bins=bins)
    print("\nProbability distribution of lytic phages:")
    for i, (bin_label, count) in enumerate(zip(bin_labels, lytic_counts)):
        print(f"  {bin_label}: {count} samples ({count / len(lytic_probs) * 100:.2f}%)")

    #
    total_counts, _ = np.histogram(probs, bins=bins)
    print("\nOverall probability distribution:")
    for i, (bin_label, count) in enumerate(zip(bin_labels, total_counts)):
        print(f"  {bin_label}: {count} samples ({count / len(probs) * 100:.2f}%)")

    
    print(f"\nsummarize the information:")
    print(f"Total sample number: {len(probs)}")
    print(f"lysogenic phage: {len(temperate_probs)} ({len(temperate_probs) / len(probs) * 100:.2f}%)")
    print(f"lytic phage: {len(lytic_probs)} ({len(lytic_probs) / len(probs) * 100:.2f}%)")

    
    plt.hist(temperate_probs, bins=50, alpha=0.7, label='lysogenic phage', color='blue')
    plt.hist(lytic_probs, bins=50, alpha=0.7, label='lytic phage', color='red')

    plt.xlabel('Prediction probability (lytic)')
    plt.ylabel('sample number')
    plt.title('Prediction probability distribution')
    plt.legend()
    plt.grid(True, alpha=0.3)

    
    stats_save_path = os.path.join(log_dir,
                                   f"{model_name}_probability_stats_{datetime.datetime.now().strftime('%m-%d-%H-%M')}.txt")
    with open(stats_save_path, 'w', encoding='utf-8') as f:
        f.write("=== Probability distribution statistics ===\n\n")
        f.write("Probability distribution of lysogenic phages:\n")
        for i, (bin_label, count) in enumerate(zip(bin_labels, temperate_counts)):
            f.write(f"  {bin_label}: {count} samples ({count / len(temperate_probs) * 100:.2f}%)\n")

        f.write("\nProbability distribution of lytic phages:\n")
        for i, (bin_label, count) in enumerate(zip(bin_labels, lytic_counts)):
            f.write(f"  {bin_label}: {count} samples ({count / len(lytic_probs) * 100:.2f}%)\n")

        f.write("\nOverall probability distribution:\n")
        for i, (bin_label, count) in enumerate(zip(bin_labels, total_counts)):
            f.write(f"  {bin_label}: {count} samples ({count / len(probs) * 100:.2f}%)\n")

        f.write(f"\nsummarize the information:\n")
        f.write(f"Total sample number: {len(probs)}\n")
        f.write(f"lysogenic phage: {len(temperate_probs)} ({len(temperate_probs) / len(probs) * 100:.2f}%)\n")
        f.write(f"lytic phage: {len(lytic_probs)} ({len(lytic_probs) / len(probs) * 100:.2f}%)\n")

    # 
    dist_save_path = os.path.join(log_dir,
                                  f"{model_name}_probability_distribution_{datetime.datetime.now().strftime('%m-%d-%H-%M')}.pdf")
    plt.savefig(dist_save_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"The probability distribution graph has been saved to: {dist_save_path}")
    print(f"The statistical results have been saved to: {stats_save_path}")

    
    stats = {
        'bins': bin_labels,
        'temperate_counts': temperate_counts,
        'lytic_counts': lytic_counts,
        'total_counts': total_counts,
        'temperate_total': len(temperate_probs),
        'lytic_total': len(lytic_probs),
        'total_samples': len(probs)
    }

    return stats


if __name__ == "__main__":
    parser = build_argparser()
    args = parser.parse_args()

    nb_classes = args._nb_classes  
    nb_samples_per_class = args._nb_samples_per_class
    seq_length = nb_classes * nb_samples_per_class

    hidden_size = args._hidden_size
    num_layers = args._num_layers

    # Create a data generator
    data_generator = PhageFeatureGenerator(
        temperate_file=args.temperate_file,
        lytic_file=args.lytic_file,
        nb_classes=nb_classes,
        nb_samples_per_class=nb_samples_per_class,
        train_ratio=args.train_ratio
    )

    # Obtain the feature dimensions from the data generator
    input_size = data_generator.feature_dim

    # create models

    model = PhageLSTMModel(
        input_size=input_size,
        hidden_size=hidden_size,
        num_layers=num_layers,
        num_classes=nb_classes,
        seq_length=seq_length,
        model_type=args.model
    ).model

    # Create a save directory
    os.makedirs(os.path.join(args.save_dir, args.model), exist_ok=True)
    os.makedirs(args.log_dir, exist_ok=True)

    if args.mode == "train":
        train(model, data_generator, args)

        # After the training is completed, load the best model for the final test
        best_model_path = os.path.join(args.save_dir, args.model, "best_model.h5")
        if os.path.exists(best_model_path):
            print("\nConduct the final test using the best model...")
            model = tf.keras.models.load_model(best_model_path)
            test(model, data_generator, args)

    elif args.mode == "test":
        # Load the pre-trained model
        model_path = os.path.join(args.save_dir, args.model, "model.h5")
        if os.path.exists(model_path):
            model = tf.keras.models.load_model(model_path)
            print(f"Loaded model from {model_path}")
        else:
            print(f"No pre-trained model found at {model_path}")
            exit(1)

        test(model, data_generator, args)