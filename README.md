# AlignmentForge

AlignmentForge is a visual software tool for biology. You use it to explore, filter, trim, and curate phylogenomic alignments.

## Access the Software

You can use AlignmentForge in two ways: as a web application or as a desktop application.

### Use the Web Application

1. Open your internet browser.
2. Go to the GitHub Pages link for this repository.
3. Use the interface to upload and analyze your alignment files.

*Note: Some features that need local file access might not operate in the web application.*

### Use the Desktop Application

The desktop application operates on Windows, macOS, and Linux. It gives you full access to all features.

1. Go to the **Releases** page of this repository.
2. Download the installer file for your operating system.
3. Install the software on your computer.
4. Open the AlignmentForge application.

## How to Build the Software

If you want to build the software from the source code, follow these steps.

### Prerequisites

You must install these tools on your computer:
* Node.js (version 20 or higher)
* Rust (latest stable version)

### Build Steps

1. Open a terminal window.
2. Clone this repository to your computer.
3. Type `npm install` and press Enter. This command installs the necessary dependencies.
4. Type `npm run dev` and press Enter. This command starts the development server.
5. Type `npm run tauri dev` and press Enter. This command starts the desktop application in development mode.

To build the final desktop application, type `npm run tauri build` and press Enter.

## Features

* **Explore Data**: See your sequence alignments in a visual grid.
* **Filter Sequences**: Remove sequences that do not meet your quality rules.
* **Trim Alignments**: Cut the ends of alignments to remove bad data.
* **Curate Data**: Select the best data for your phylogenetic analysis.

## License

Please refer to the LICENSE file in this repository.
